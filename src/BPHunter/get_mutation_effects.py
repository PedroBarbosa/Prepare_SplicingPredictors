import sys
import pandas as pd 
import numpy as np

df = pd.read_csv(sys.argv[1], sep="\t")
df['CHROM'] = df.CHROM.str.replace("chr", "")
df.columns = [x.strip() for x in df.columns]

chroms = list(set([x for x in df.CHROM.unique()]))
vcf_header = "##fileformat=VCFv4.2\n"
vcf_header += "{}\n".format('\n'.join(["##contig=<ID=" + x + ">" for x in sorted(chroms)]))
vcf_header += "##INFO=<ID=BPHunter,Number=1,Type=String,Description=\"BPHunter score, genome-wide detection of human intronic variants that disrupt the branchpoint signal.\">\n"
vcf_header += "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO"
print(vcf_header)
for i, row in df.iterrows():
    chrom = row[0]
    pos = row[1]
    v_id = row[2]
    ref = row[3]
    alt = row[4]
    info = "BPHunter={};".format(row[-1])
    print('{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}'.format(chrom, pos, v_id, ref, alt, '.', '.', info))