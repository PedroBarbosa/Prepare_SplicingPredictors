We use the ESRseq hexamer scores retrieved from the [original publication](https://genome.cshlp.org/content/21/8/1360), supplementary table 1. The scores were saved in `scores.tsv` file.

The variant effect is calculated as in [here](https://onlinelibrary.wiley.com/doi/10.1002/humu.22428]). Briefly, the hexamers overlapping the position of the wildtype and mutated sequences are retrieved. Then, for each group the individual hexamer scores are summed to calculated the total ESRseq score per group. The change is calculated by subtracting the ESRseq mutated - ESRseq WT (in the paper it was a different orientation).

- Generate input sequences from VCF:\
`vcf2seq infile.vcf.gz hg19_genome.fa out --esrseq`

- Calculate mutation effect and generate new VCF (`scores.tsv` should exist in script directory):\
`python get_mutation_effect out_ESRseq.fa | bcftools sort | bgzip > ESRseq_to_annotate.vcf.gz`\

- Optionally, annotate the `ESRseq` scores into the main VCF file (to benchmark with other tools) using `vcfanno` or `bcftools annotate`
