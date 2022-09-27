Step to generate DSSP finder results:

- Generate input sequences from VCF:\
`vcf2seq infile_donor.vcf.gz hg19_genome.fa out -ss donor --dssp`\
`vcf2seq infile_acceptor.vcf.gz hg19_genome.fa out -ss acceptor --dssp`

- Clone the [official repo](https://github.com/DSSP-github/DSSP), make sure the dependencies are installed, and run the donor and acceptor models:\
`python DS_DSSP.py -I out_DSSP_donor.fa -O dssp_donor_out.tsv`
`python AS_DSSP.py -I out_DSSP_acceptor.fa -O dssp_acceptor_out.tsv`\

- To get a ready-to-annotate VCF file, just run:\
`python get_mutation_effect.py dssp_donor_out.tsv| bcftools sort | bgzip > dssp_donor_to_annotate.vcf.gz`\
`python get_mutation_effect.py dssp_acceptor_out.tsv| bcftools sort | bgzip > dssp_acceptor_to_annotate.vcf.gz`

- Optionally, annotate the `DSSP` scores into the main VCF file (to benchmark with other tools) using `vcfanno` or `bcftools annotate`
