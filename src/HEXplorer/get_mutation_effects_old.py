import sys
import pandas as pd 
import numpy as np

df = pd.read_csv(sys.argv[1], index_col=False, sep=";", low_memory=False).rename_axis('idx').reset_index()
df.columns = [x.strip() for x in df.columns]

mut_map = pd.read_csv(sys.argv[2], sep="\t", names=['v_id', 'mut_idx'])
mut_map['hexamer_idx'] = mut_map.mut_idx.apply(lambda x: range(x-5, x+6))
mut_map = mut_map.explode('hexamer_idx')

df = pd.merge(mut_map, df, left_on='hexamer_idx', right_on='idx', how='left')
df.drop(['idx'], axis=1, inplace=True)

def get_diff(group: pd.DataFrame):
    max_diff = group.loc[group.apply(lambda x: x.HZEI_mt - x.HZEI_wt, axis=1).abs().idxmax()]
    return str(max_diff.mut_idx - max_diff.hexamer_idx) + "|" + str(round(max_diff.HZEI_mt - max_diff.HZEI_wt, 3))

df[['HZEI_wt', 'HZEI_mt']] = df[['HZEI_wt', 'HZEI_mt']].astype(float)
df = df.groupby(['v_id']).apply(get_diff).reset_index()
df = df.rename(columns={0: "pred"})

chroms = list(set([x.split("_")[0] for x in df.v_id.unique()]))
vcf_header = "##fileformat=VCFv4.2\n"
vcf_header += "{}\n".format('\n'.join(["##contig=<ID=" + x + ">" for x in sorted(chroms)]))
vcf_header += "##INFO=<ID=HEXplorer,Number=1,Type=String,Description=\"Max HZEI difference between mutated and WT in the hexamers overlapping the variant. Format: Distance_to_variant|Prediction_difference.\">\n"
vcf_header += "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO"
print(vcf_header)
for i, row in df.iterrows():
    fields = row.v_id.split("_")
    chrom = fields[0]
    pos = fields[1]
    ref = fields[2]
    alt = fields[3]
    info = "HEXplorer={};".format(row.pred)
    print('{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}'.format(chrom, pos, row.v_id, ref, alt, '.', '.', info))