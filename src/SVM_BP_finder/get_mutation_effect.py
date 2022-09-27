import sys
from collections import defaultdict

preds = defaultdict(dict) 
for _line in sys.stdin:

    if isinstance(_line, str):
        
        line = _line.rstrip()
        v_id = '_'.join(line.split()[0].split("_")[:-1])
        group = line.split()[0].split("_")[-1]
        preds[v_id].update({group: float(line.split()[-1])})

score_diff = {}
for v_id, scores in preds.items():
    if len(scores) != 2:
        raise ValueError('Variant {} does not have a score for WT and mutated seq'.format(v_id))

    score_diff[v_id] = round((scores['Mutated'] - scores['WT']) / scores['WT'], 3)

chroms = list(set([x.split("_")[0] for x in score_diff.keys()]))
vcf_header = "##fileformat=VCFv4.2\n"
vcf_header += "{}\n".format('\n'.join(["##contig=<ID=" + x + ">" for x in sorted(chroms)]))
vcf_header += "##INFO=<ID=SVM_BP_finder,Number=1,Type=Float,Description=\"Difference between the score for the strongest branchpoint position at the reference and mutated sequence\">\n"
vcf_header += "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO"
print(vcf_header)
for v_id, value in score_diff.items():
    fields = v_id.split("_")
    chrom = fields[0]
    pos = fields[1]
    ref = fields[2]
    alt = fields[3]
    info = "SVM_BP_finder={}".format(str(value))
    print('{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}'.format(chrom, pos, v_id, ref, alt, '.', '.', info))