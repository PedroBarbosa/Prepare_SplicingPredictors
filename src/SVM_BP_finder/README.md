Step to generate SVM-BP finder results:

- Generate input sequences from VCF:\
`vcf2seq infile.vcf.gz hg19_genome.fa out --svm_bp_finder`

- Run the model by uploading the `out_SVM_BP_finder.fa` file to [available server](http://regulatorygenomics.upf.edu/Software/SVM_BP/) and saving the results into a new file (e.g. svm\_bp\_out.tsv)

- Select the best BP per sequence, as suggested by the authors, and automatically generate a VCF :\
`perl calculate_best_BP_per_intron.pl < svm_bp_out.tsv | sort -k1 | python get_mutation_effect.py | bcftools sort | bgzip > SVM_BP_finder.vcf.gz`
    
    - It is crucial to sort the `calculate_best_BP_per_intron` output so that the outputs from the same variant (reference and mutated sequence) come 
    together. 

- Optionally, annotate the `SVM_BP_finder` scores into the main VCF file (to benchmark with other tools) using `vcfanno` or `bcftools annotate`
