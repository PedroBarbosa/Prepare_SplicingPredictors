import sys
import pandas as pd 
import numpy as np

df = pd.read_csv(sys.argv[1], index_col=False, sep=";", decimal=',', low_memory=False)
df.columns = [x.strip() for x in df.columns]

mut_map = pd.read_csv(sys.argv[2], sep="\t", names=['v_id', 'mut_idx'])
mut_map['hexamer_idx'] = mut_map.mut_idx.apply(lambda x: range(x-5, x+6))
mut_map = mut_map.explode('hexamer_idx')
v_ids = mut_map.v_id.unique()

df_ref = df[df.Sequence == "reference"]
df_alt = df[df.Sequence == "alternative"]

hexamer = 11
v_idx = 0
out = []
    
assert len(df_ref) == len(df_alt)

complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
for i in range(0, len(df_ref), hexamer):

    v_id = v_ids[v_idx]

    ref_allele = v_id.split("_")[2]
    alt_allele = v_id.split("_")[3]
    strand = v_id.split("_")[-1]
    if strand == "rev":
        ref_allele = complement[ref_allele]
        alt_allele = complement[alt_allele]
        
    batch_ref = df_ref.iloc[i:i+hexamer].head(6)
    batch_alt = df_alt.iloc[i:i+hexamer].head(6)

    assert batch_ref.iloc[0].seq9[5] == ref_allele
    assert batch_alt.iloc[0].seq9[5] == alt_allele
    
    df = pd.merge(batch_ref[['durchzahl', 'endhex']] , batch_alt[['durchzahl', 'endhex']], on='durchzahl', how='left')
    max_idx = df.apply(lambda x: x.endhex_y - x.endhex_x, axis=1).abs().idxmax()
    max_diff = round(df.loc[max_idx].endhex_y - df.loc[max_idx].endhex_x, 3)
    hexamer_pos = max_idx + 1
    out.append([v_id, hexamer_pos, max_diff])
    v_idx += 1
    
df = pd.DataFrame(out, columns=['v_id', 'hexamer_pos', 'max_diff'])


chroms = list(set([x.split("_")[0] for x in df.v_id.unique()]))
vcf_header = "##fileformat=VCFv4.2\n"
vcf_header += "{}\n".format('\n'.join(["##contig=<ID=" + x + ">" for x in sorted(chroms)]))
vcf_header += "##INFO=<ID=HEXplorer,Number=1,Type=String,Description=\"Max HZEI difference between mutated and WT in the hexamers overlapping the variant. Format: Position_in_hexamer|Prediction_difference.\">\n"
vcf_header += "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO"
print(vcf_header)
for i, row in df.iterrows():
    fields = row.v_id.split("_")
    chrom = fields[0]
    pos = fields[1]
    ref = fields[2]
    alt = fields[3]
    info = "HEXplorer={};".format(str(row.hexamer_pos) + "|" + str(row.max_diff))
    print('{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}'.format(chrom, pos, row.v_id, ref, alt, '.', '.', info))