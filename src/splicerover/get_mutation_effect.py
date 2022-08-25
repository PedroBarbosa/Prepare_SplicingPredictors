import sys
from Bio import SeqIO
import pandas as pd 

# Fasta input
handle = open(sys.argv[1])
map_fasta_ids = {}
for i, record in enumerate(SeqIO.parse(handle, "fasta")):
    v_id = '_'.join(record.id.split()[0].split("_"))
    map_fasta_ids[i] = v_id

map_fasta_ids = pd.DataFrame.from_dict(map_fasta_ids, orient='index', 
                                       columns=['seq_id']).rename_axis('id').reset_index()

# Preds
preds_df = pd.read_csv(sys.argv[2], index_col=False, names=['id','pos', 'seq_motif', 'pred'], sep=",")
preds_df['pos'] += 1
df = pd.merge(preds_df, map_fasta_ids, on='id', how='left')

df[['v_id', 'group']] = df.seq_id.str.rsplit('_', n=1, expand=True)
df.drop(columns='seq_id', inplace=True)

def get_diff(group: pd.DataFrame):

    nrows = group.shape[0]
    if nrows == 1:
        
        if group.iloc[0,:].group == "Mutated":
            return group.iloc[0,:].pred.item()
        
        if group.iloc[0,:].group == "WT":
            return - group.iloc[0,:].pred.item()
    
    elif nrows == 2:    

        mut = group[group.group == "Mutated"].pred.item()
        wt = group[group.group == "WT"].pred.item()
        return mut - wt

    else:
        raise ValueError('Too many rows in groupby')

df['pred'] = df.pred.astype(float)
df = df.groupby(['v_id', 'pos']).apply(get_diff).reset_index()
df = df.rename(columns={0: 'diff'})
df = df.loc[df.groupby('v_id')['diff'].apply(lambda x: x.abs().idxmax())]

chroms = list(set([x.split("_")[0] for x in df.v_id.unique()]))
vcf_header = "##fileformat=VCFv4.2\n"
vcf_header += "{}\n".format('\n'.join(["##contig=<ID=" + x + ">" for x in sorted(chroms)]))
vcf_header += "##INFO=<ID=SpliceRover,Number=1,Type=String,Description=\"Max difference at a given position between the predictions for the reference and mutated sequences. Format: Position|Prediction_difference.\">\n"
vcf_header += "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO"
print(vcf_header)
for i, row in df.iterrows():
    fields = row.v_id.split("_")
    chrom = fields[0]
    pos = fields[1]
    ref = fields[2]
    alt = fields[3]
    info = "SpliceRover={}|{};".format(str(row.pos), str(round(row['diff'], 3)))
    print('{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}'.format(chrom, pos, row.v_id, ref, alt, '.', '.', info))