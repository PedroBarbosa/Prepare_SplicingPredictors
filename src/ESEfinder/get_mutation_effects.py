import sys
import pandas as pd 
import numpy as np

# ESEfinder output
file = open(sys.argv[1], 'r')

target_lines = False
df = []
for line in file:
    if line.startswith('ESEfinder result:'):
        target_lines = True
        continue
    if target_lines and len(line.split("\t")) > 1:
        df.append(line.rstrip().split("\t"))

df = pd.DataFrame.from_records(df[1:], columns=df[0])
df[['v_id', 'group']] = df.Seq.str.rsplit('_', n=1, expand=True)
df.drop(columns='Seq', inplace=True)

def get_diff(group: pd.DataFrame):
    assert group.shape[0] == 2, "There should be WT and Mut seqs for the given v_id, position and motif"
    wt = group[group.group == "WT"].Score.item()
    mut = group[group.group == "Mutated"].Score.item()
    return mut - wt

def select_highest(group: pd.DataFrame):
    max_diff = group.loc[group['diff'].abs().idxmax()]
    return str(max_diff.Motif) + "|" + str(round(max_diff['diff'], 3))

df['Score'] = df.Score.astype(float)
diffs = df.groupby(['v_id', 'Motif', 'Position']).apply(get_diff).reset_index()
diffs = diffs.rename(columns={0: "diff"})
diffs['Motif'] = diffs.Motif.str.replace(" ", "")
results = diffs.groupby('v_id').apply(select_highest).reset_index()
results = results.rename(columns={0: "pred"})

chroms = list(set([x.split("_")[0] for x in df.v_id.unique()]))
vcf_header = "##fileformat=VCFv4.2\n"
vcf_header += "{}\n".format('\n'.join(["##contig=<ID=" + x + ">" for x in sorted(chroms)]))
vcf_header += "##INFO=<ID=ESEfinder,Number=1,Type=String,Description=\"Max ESE score difference between mutated and WT for any given SR weight matrix. Format: Motif_with_highest_difference|Score_difference.\">\n"
vcf_header += "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO"
print(vcf_header)
for i, row in results.iterrows():
    fields = row.v_id.split("_")
    chrom = fields[0]
    pos = fields[1]
    ref = fields[2]
    alt = fields[3]
    info = "ESEfinder={};".format(row.pred)
    print('{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}'.format(chrom, pos, row.v_id, ref, alt, '.', '.', info))