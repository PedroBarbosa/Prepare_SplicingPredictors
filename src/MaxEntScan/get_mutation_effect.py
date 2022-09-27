import sys
from Bio import SeqIO
import pandas as pd
from maxentpy import maxent
from maxentpy.maxent import load_matrix5, load_matrix3

fasta = sys.argv[1]
ss = sys.argv[2]

assert ss in ['donor', 'acceptor'], 'Second argument must tell what kind of splice site to score'

# Fasta input
handle = open(fasta)
to_df = []
for i, record in enumerate(SeqIO.parse(handle, "fasta")):
    v_id = '_'.join(record.id.split()[0].split("_")[:-1])
    group = record.id.split()[0].split("_")[-1]
    to_df.append([v_id, group, str(record.seq)])

df = pd.DataFrame.from_records(to_df, columns=['v_id', 'group', 'seq'])

def score(wt: str, mut:str, ss:str, matrix, i: int = None):
    
    if ss == "donor":
        wt =  wt[:3].lower() + wt[3:].upper()
        mut = mut[:3].lower() + mut[3:].upper()
        score_wt = maxent.score5(wt, matrix=matrix)
        score_mut = maxent.score5(mut, matrix=matrix)
        #dist_to_mut = 5 - i 
        #mut_position = 9 - i
        
    else:
        wt =  wt[:18].lower() + wt[18:20].upper() + wt[20:].lower()
        mut = mut[:18].lower() + mut[18:20].upper() + mut[20:].lower()

        score_wt = maxent.score3(wt, matrix=matrix)
        score_mut = maxent.score3(mut, matrix=matrix)
        #dist_to_mut = 3 - i 
        #mut_position = 23 - i
    
    return score_wt, score_mut

def get_maxentscan_result(group: pd.DataFrame, matrix, ss):
    wt = group[group.group == "WT"].seq.item()
    mut = group[group.group == "Mutated"].seq.item()
   
    subseq_size = 9 if ss == "donor" else 23
    subseqs_wt = [wt[i:i+subseq_size] for i in range(len(wt)-subseq_size + 1)]
    subseqs_mut = [mut[i:i+subseq_size] for i in range(len(mut)-subseq_size + 1)]
    
    # Fixed pos (new ss, loss ss)
    _idx = 4 if ss == "donor" else 3
    seq_wt = subseqs_wt[_idx: _idx + 2]
    seq_mut = subseqs_mut[_idx: _idx + 2] 
    out = []
    for wt, mut in zip(seq_wt, seq_mut):
        score_wt, score_mut = score(wt, mut, ss, matrix)
        out.append(score_mut - score_wt)
    fixed_diff = max(out, key=abs)
 
    # Sliding approach
    out = []
    for i in range(len(subseqs_wt)):
        wt = subseqs_wt[i]
        mut = subseqs_mut[i]
        score_wt, score_mut = score(wt, mut, ss, matrix, i=i)
        
        out.append([score_wt, score_mut, wt, mut])
        #out.append([score_mut, score_wt, dist_to_mut, mut_position, wt, mut])
    
    out = pd.DataFrame.from_records(out, columns=['score_wt', 'score_mut', 'wt', 'mut'])
    #out = pd.DataFrame.from_records(out, columns=['score_mut', 'score_wt', 'dist_to_mutation', 'mut_position', 'wt', 'mut'])

    sliding_diff = max(out.score_mut) - max(out.score_wt)
    return str(round(max([fixed_diff, sliding_diff], key=abs), 3))

 
matrix = load_matrix5() if ss == "donor" else load_matrix3()
results = df.groupby('v_id').apply(get_maxentscan_result, matrix=matrix, ss=ss)

chroms = list(set([x.split("_")[0] for x in results.index]))
vcf_header = "##fileformat=VCFv4.2\n"
vcf_header += "{}\n".format('\n'.join(["##contig=<ID=" + x + ">" for x in sorted(chroms)]))
vcf_header += "##INFO=<ID=MaxEntScan,Number=1,Type=String,Description=\"Max difference through a sliding window between the predictions for the reference and mutated sequences. Format: Distance_to_mutation|Prediction_difference.\">\n"
vcf_header += "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO"
print(vcf_header)
for v_id, pred in results.items():

    fields = v_id.split("_")
    chrom = fields[0]
    pos = fields[1]
    ref = fields[2]
    alt = fields[3]
    info = "MaxEntScan={};".format(pred)
    print('{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}'.format(chrom, pos, v_id, ref, alt, '.', '.', info))

