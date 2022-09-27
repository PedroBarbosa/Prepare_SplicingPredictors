import sys
import pandas as pd 
import os 
import numpy as np
from collections import defaultdict
from Bio import SeqIO
from pyparsing import java_style_comment

# ESRseq sequences
file = open(sys.argv[1], 'r')

scores = os.path.join(os.path.dirname(os.path.realpath(__file__)), 'scores.tsv')
hexamers = pd.read_csv(scores, sep="\t", names=['hexamer', 'score'])

results =  defaultdict(dict)
for i, record in enumerate(SeqIO.parse(file, "fasta")):
    v_id = '_'.join(record.id.split()[0].split("_")[:-1])
    group = record.id.split()[0].split("_")[-1]
    
    k_mer_size = 6
    seq = str(record.seq)
    total_esrseq_score = 0
    for i in range(6):

        kmer = seq[i:i+k_mer_size]
        _score = hexamers[hexamers.hexamer == kmer]
        if not _score.empty:
            total_esrseq_score += _score.score.item()
    
    results[v_id].update({group: total_esrseq_score})

score_diff = {}
for v_id, _preds in results.items():
    if len(_preds) != 2:
        raise ValueError('Variant {} does not have a score for WT and mutated seq'.format(v_id))
    score_diff[v_id] = round(_preds['Mutated'] - _preds['WT'], 3)

chroms = list(set([x.split("_")[0] for x in score_diff.keys()]))
vcf_header = "##fileformat=VCFv4.2\n"
vcf_header += "{}\n".format('\n'.join(["##contig=<ID=" + x + ">" for x in sorted(chroms)]))
vcf_header += "##INFO=<ID=ESRseq,Number=1,Type=Float,Description=\"Delta ESRseq scores: total ESRseq scores for the mutated sequence - total ESRseq scores for the wildtype sequence\">\n"
vcf_header += "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO"
print(vcf_header)
for v_id, value in score_diff.items():
    fields = v_id.split("_")
    chrom = fields[0]
    pos = fields[1]
    ref = fields[2]
    alt = fields[3]
    info = "ESRseq={}".format(str(value))
    print('{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}'.format(chrom, pos, v_id, ref, alt, '.', '.', info))