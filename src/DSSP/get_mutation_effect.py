import sys
from collections import defaultdict

preds = defaultdict(dict)
dssp_preds = open(sys.argv[1])
for line in dssp_preds:
    v_id = '_'.join(line[1:].split()[0].split("_")[:-1])
    group = line.split()[0].split("_")[-1]
    
    score = next(dssp_preds)
    preds[v_id].update({group: float(score)})

score_diff = {}
for v_id, scores in preds.items():
    if len(scores) != 2:
        raise ValueError('Variant {} does not have a score for WT and mutated seq'.format(v_id))

    score_diff[v_id] = round(scores['Mutated'] - scores['WT'], 3)
    
chroms = list(set([x.split("_")[0] for x in score_diff.keys()]))
vcf_header = "##fileformat=VCFv4.2\n"
vcf_header += "{}\n".format('\n'.join(["##contig=<ID=" + x + ">" for x in sorted(chroms)]))
vcf_header += "##INFO=<ID=DSSP,Number=1,Type=Float,Description=\"Score difference between mutated and wiltype sequence. Donor pos (71/72); Acceptor pos (69/70)\">\n"
vcf_header += "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO"
print(vcf_header)
for v_id, value in score_diff.items():
    fields = v_id.split("_")
    chrom = fields[0]
    pos = fields[1]
    ref = fields[2]
    alt = fields[3]
    info = "DSSP={}".format(str(value))
    print('{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}'.format(chrom, pos, v_id, ref, alt, '.', '.', info))