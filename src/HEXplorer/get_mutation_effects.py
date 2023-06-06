import sys
import pandas as pd 
import numpy as np

df = pd.read_csv(sys.argv[1], index_col=False, sep=";", decimal=",", low_memory=False)
df.columns = [x.strip() for x in df.columns]

mut_map = pd.read_csv(sys.argv[2], sep="\t", names=['v_id', 'mut_idx'])
mut_map['hexamer_idx'] = mut_map.mut_idx.apply(lambda x: range(x-5, x+6))
mut_map = mut_map.explode('hexamer_idx')
v_ids = mut_map.v_id.unique()

# Because HEXplorer doesn't report scores for the first five nucleotides, we adjust the mutation indices accordingly
mut_map['mut_idx'] = mut_map.mut_idx - 5 
complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}

out = []
for v_id, group in mut_map.groupby('v_id', sort=False):
  
    mut_index = group.iloc[0].mut_idx

    _df = df.loc[mut_index-5:mut_index+6].copy()
    variant_nuc = _df.iloc[5].indexNT
    ref_allele = v_id.split("_")[2]
    alt_allele = v_id.split("_")[3]
    strand = v_id.split("_")[-1]
    if strand == "rev":
        ref_allele = complement[ref_allele]
        alt_allele = complement[alt_allele]

    assert variant_nuc == ref_allele, f"Something went wrong for variant {v_id}. indexNT doesn't match the reference allele."

    _df['HZEI_diff'] = _df['mt_HZEI'] - _df['wt_HZEI']
    final_score = _df['HZEI_diff'].sum()
    out.append([v_id, final_score])

df = pd.DataFrame(out, columns=['v_id', 'pred'])
chroms = list(set([x.split("_")[0] for x in v_ids]))
vcf_header = "##fileformat=VCFv4.2\n"
vcf_header += "{}\n".format('\n'.join(["##contig=<ID=" + x + ">" for x in sorted(chroms)]))
vcf_header += "##INFO=<ID=HEXplorer,Number=1,Type=String,Description=\"Sum of the HZEI differences between mutated and WT in the hexamers overlapping the variant.\">\n"
vcf_header += "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO"
print(vcf_header)
for i, row in df.iterrows():
    fields = row.v_id.split("_")
    chrom = fields[0]
    pos = fields[1]
    ref = fields[2]
    alt = fields[3]
    info = "HEXplorer={};".format(round(row.pred,3))
    print('{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}'.format(chrom, pos, row.v_id, ref, alt, '.', '.', info))