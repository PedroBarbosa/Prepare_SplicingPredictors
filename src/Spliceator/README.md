Step to generate Spliceator finder results:

- Generate input sequences from VCF:\
`vcf2seq infile_donors.vcf.gz hg19_genome.fa out --spliceator --ss donor`\
`vcf2seq infile_acceptors.vcf.gz hg19_genome.fa out --spliceator --ss acceptor`

- Clone the [official repo](https://git.unistra.fr/nscalzitti/spliceator), make sure the dependecnies are installed, and run the model by generating predictions for all positions in the input sequence:\
`./Spliceator --all True -f Spliceator_donor.fa -o spliceator_donor_out.csv'`\
`./Spliceator --all True -f Spliceator_acceptor.fa -o spliceator_acceptor_out.csv'`

- The model is slow. The output is multi-column file separated by ";", with predictions for all input positions. To get a ready-to-annotate VCF file, just run:\
`python get_mutation_effect.py spliceator_donor_out.csv Donor | bcftools sort | bgzip > spliceator_donor_to_annotate.vcf.gz`\
`python get_mutation_effect.py spliceator_acceptor_out.csv Acceptor | bcftools sort | bgzip > spliceator_acceptor_to_annotate.vcf.gz`

- Optionally, annotate the `Spliceator` scores into the main VCF file (to benchmark with other tools) using `vcfanno` or `bcftools annotate`
