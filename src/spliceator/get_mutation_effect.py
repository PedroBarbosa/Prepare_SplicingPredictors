import sys
import pandas as pd 

df = pd.read_csv(sys.argv[1], sep=";", low_memory=False)
splice_site = sys.argv[2]

assert splice_site in ['Donor', 'Acceptor'], "Second argument must be the type of splice site"

df = df[(df.SS_type == splice_site) & (df.ID != "ID")]
df[['v_id', 'group']] = df.ID.str.rsplit('_', n=1, expand=True)
df.drop(columns='ID', inplace=True)

def get_diff(group: pd.DataFrame):
    assert group.shape[0] == 2
    mut = group[group.group == "Mutated"].Score.item()
    wt = group[group.group == "WT"].Score.item()
    return mut - wt

df['Score'] = df.Score.astype(float)
df = df.groupby(['v_id', '#']).apply(get_diff).reset_index()
df = df.rename(columns={0: "diff"})
df = df.loc[df.groupby('v_id')['diff'].apply(lambda x: x.abs().idxmax())]

chroms = list(set([x.split("_")[0] for x in df.v_id.unique()]))
vcf_header = "##fileformat=VCFv4.2\n"
vcf_header += "{}\n".format('\n'.join(["##contig=<ID=" + x + ">" for x in sorted(chroms)]))
vcf_header += "##INFO=<ID=Spliceator,Number=1,Type=String,Description=\"Max difference at a given position between the predictions for the reference and mutated sequences. Format: Position|Prediction_difference.\">\n"
vcf_header += "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO"
print(vcf_header)
for i, row in df.iterrows():
    fields = row.v_id.split("_")
    chrom = fields[0]
    pos = fields[1]
    ref = fields[2]
    alt = fields[3]
    info = "Spliceator={}|{};".format(row['#'].replace("#", ""), str(round(row['diff'], 3)))
    print('{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}'.format(chrom, pos, row.v_id, ref, alt, '.', '.', info))