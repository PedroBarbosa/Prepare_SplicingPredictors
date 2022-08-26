Step to generate MaxEntScan results:

- Generate input sequences (to be tested in a sliding window) from VCF:\
`vcf2seq infile_donor.vcf.gz hg19_genome.fa out -ss donor --maxentscan`\
`vcf2seq infile_acceptor.vcf.gz hg19_genome.fa out --ss acceptor --maxentscan`

- Run maxentscan (make sure the package maxentscanpy is installed) and save the max effect of the mutation:\
`python get_mutation_effect.py new_donor_MaxEntScan_donor.fa donor | bcftools sort | bgzip > maxentscan_donor_to_annotate.vcf.gz`\
`python get_mutation_effect.py new_acceptor_MaxEntScan_acceptor.fa acceptor| bcftools sort | bgzip > maxentscan_acceptor_to_annotate.vcf.gz`

- Optionally, annotate the `MaxEntScan` scores into the main VCF file (to benchmark with other tools) using `vcfanno` or `bcftools annotate`