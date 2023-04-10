import sys
import pandas as pd 

if sys.argv[2]:
    assert sys.argv[2] == "shift", "If providing an additional argument, it must be 'shift', meaning that variants with score 0 will be shifted to 1"
    shift = True
else:
    shift = False
   
df = pd.read_csv(sys.argv[1], sep="\t")
df['CHROM'] = df.CHROM.str.replace("chr", "")
df.columns = [x.strip() for x in df.columns]

df = df.drop(columns='BPHUNTER_SCORE')
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
    score = row[-1]
    if shift and score == 0:
        score = 1
        
    info = "BPHunter={};".format(score)
    print('{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}'.format(chrom, pos, v_id, ref, alt, '.', '.', info))