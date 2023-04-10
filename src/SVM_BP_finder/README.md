Step to generate SVM-BP finder results:

- Generate input sequences from VCF (500bp upstream of acceptor site):\
`vcf2seq infile.vcf.gz hg19_genome.fa out --svm_bp_finder`

- Download the [official repo](https://github.com/comprna/SVM-BPfinder), make sure the dependencies are installed (adapt prints for Python3, if necessary), and run the model accordingly:\
`./svm_bpfinder.py -i input_SVM_BP_finder.fa -s Hsap -l 500 -d 50 > svm_bp_out.tsv`

- Select the best BP per sequence, as suggested by the authors, and automatically generate a VCF :\
`perl calculate_best_BP_per_intron.pl < svm_bp_out.tsv | sort -k1,1 | python get_mutation_effect.py | bcftools sort | bgzip > SVM_BP_finder.vcf.gz`
    
    - It is crucial to sort the `calculate_best_BP_per_intron` output so that the outputs from the same variant (reference and mutated sequence) come 
    together. 

- Optionally, annotate the `SVM_BP_finder` scores into the main VCF file (to benchmark with other tools) using `vcfanno` or `bcftools annotate`
