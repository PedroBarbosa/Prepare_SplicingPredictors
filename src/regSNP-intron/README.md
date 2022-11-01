Step to generate regSNP-intron results:

- Select only SNVs and generate proper input format:\
`bcftools view -v snps file.vcf | grep -v "^#" | cut -f1,2,4,5 > regSNP_intron_input.tsv`

- If vcf does not have 'chr' in chromosome names, please add that.\
`bcftools view -v snps file.vcf | grep -v "^#" | sed 's/^/chr/g' | cut -f1,2,4,5 > regSNP_intron_input.tsv`

- Go to the [official web app](https://regsnps-intron.ccbb.iupui.edu/), and run the model. When finished, download the predictions `snp.prediction.txt`.

- To get a ready-to-annotate VCF file, just run:\
`python get_mutation_effect.py snp.prediction.txt | bcftools sort | bgzip > regSNP_intron_to_annotate.vcf.gz`

- Optionally, annotate the `regSNP_intron` scores into the main VCF file (to benchmark with other tools) using `vcfanno` or `bcftools annotate`
