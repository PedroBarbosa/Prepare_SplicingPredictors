import sys
from Bio import SeqIO
from collections import defaultdict

# Fasta input
handle = open(sys.argv[1])

# Preds
preds = [x for x in open(sys.argv[2], 'r')]

scores =  defaultdict(dict)
for i, record in enumerate(SeqIO.parse(handle, "fasta")):
    v_id = '_'.join(record.id.split()[0].split("_")[:-1])
    group = record.id.split()[0].split("_")[-1]
    
    scores[v_id].update({group: float(preds[i].rstrip())})
 
score_diff = {}
for v_id, _preds in scores.items():
    if len(_preds) != 2:
        raise ValueError('Variant {} does not have a score for WT and mutated seq'.format(v_id))
    score_diff[v_id] = _preds['Mutated'] - _preds['WT']

chroms = list(set([x.split("_")[0] for x in score_diff.keys()]))
vcf_header = "##fileformat=VCFv4.2\n"
vcf_header += "{}\n".format('\n'.join(["##contig=<ID=" + x + ">" for x in sorted(chroms)]))
vcf_header += "##INFO=<ID=Splice2Deep,Number=1,Type=Float,Description=\"Difference between the class prediction at the mutated and reference sequence. 0=no diff. -1=splice site lost with mutation. -1=splice site gained with mutation\">\n"
vcf_header += "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO"
print(vcf_header)
for v_id, value in score_diff.items():
    fields = v_id.split("_")
    chrom = fields[0]
    pos = fields[1]
    ref = fields[2]
    alt = fields[3]
    info = "Splice2Deep={}".format(str(value))
    print('{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}'.format(chrom, pos, v_id, ref, alt, '.', '.', info))