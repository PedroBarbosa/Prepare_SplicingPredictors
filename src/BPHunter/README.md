Step to generate BPHunter results:

- No need to generate input sequences, as the software accepts VCF format. However, VCF should come without the header.\
`zcat file.vcf.gz | grep -v "^#" | cut -f1,2,3,4,5 > BPhunter_input.tsv`

- Go to the [official webapp](https://hgidsoft.rockefeller.edu/BPHunter/index.php), and run the model accordingly.

- Download the output and run our script to generate a VCF :\
`python get_mutation_effect.py BPHunter_out.tsv | bcftools sort | bgzip > BPHunter.vcf.gz`

- Optionally, annotate the `BPHunter` scores into the main VCF file (to benchmark with other tools) using `vcfanno` or `bcftools annotate`.
