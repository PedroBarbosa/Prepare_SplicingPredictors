import sys
import pandas as pd

# Raw predictions
preds = pd.read_csv(sys.argv[1], sep="\t", usecols=['#chrom', 'pos', 'ref', 'alt', 'disease', 'prob'])

chroms = list(set([x.replace('chr', '') for x in preds['#chrom'].unique()]))

vcf_header = "##fileformat=VCFv4.2\n"
vcf_header += "{}\n".format('\n'.join(["##contig=<ID=" + x + ">" for x in sorted(chroms)]))
vcf_header += "##INFO=<ID=regSNP_intron,Number=1,Type=String,Description=\"Output of regSNP_intron model. Format: Disease_class|Damage_probability.\">\n"
vcf_header += "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO"
print(vcf_header)

for _, value in preds.iterrows():

    chrom = value['#chrom']
    pos = str(value.pos)
    ref = value.ref
    alt = value.alt
    info = "regSNP_intron={}|{}".format(value.disease, str(round(value.prob, 5)))
    print('{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}'.format(chrom, pos, '{}'.format(chrom + "_" + pos + "_" + ref + "_" + alt), ref, alt, '.', '.', info))