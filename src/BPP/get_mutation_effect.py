import sys
from collections import defaultdict

file = open(sys.argv[1], 'r')
preds = defaultdict(dict) 

for _line in file:
  
    # header
    #id,bps,bp_pos,sc_bps,sc_ppt,sc,zsc_bps,zsc_ppt,zsc
    # BPS score (sc_bps) and z-score BPS (zsc_bps) are being selected
    if not _line.startswith('#'):
        
        line = _line.rstrip()
  
        v_id = '_'.join(line.split()[0].split("_")[:-1])
        v_id = v_id.replace('>', '')
        group = line.split()[0].split("_")[-1]
        
        # BPS score
        #preds[v_id].update({group: float(line.split()[3])})
        
        # BPP overall score
        preds[v_id].update({group: float(line.split()[5])})
        
        # BPS z-score 
        #preds[v_id].update({group: float(line.split()[6])})
        
        # BPP overall z-score 
        #preds[v_id].update({group: float(line.split()[8])})


score_diff = {}
for v_id, scores in preds.items():
    if len(scores) != 2:
        raise ValueError('Variant {} does not have a score for WT and mutated seq'.format(v_id))
 
    score_diff[v_id] = (scores['Mutated'] - scores['WT']) / scores['WT']

chroms = list(set([x.split("_")[0] for x in score_diff.keys()]))
vcf_header = "##fileformat=VCFv4.2\n"
vcf_header += "{}\n".format('\n'.join(["##contig=<ID=" + x + ">" for x in sorted(chroms)]))
vcf_header += "##INFO=<ID=BPP,Number=1,Type=Float,Description=\"Delta score between strongest branchpoint position at the reference and mutated sequence ((score_mut - score_wt) / score_wt).\">\n"
vcf_header += "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO"
print(vcf_header)
for v_id, value in score_diff.items():
    fields = v_id.split("_")
    chrom = fields[0]
    pos = fields[1]
    ref = fields[2]
    alt = fields[3]
    info = "BPP={}".format(str(value))
    print('{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}'.format(chrom, pos, v_id, ref, alt, '.', '.', info))