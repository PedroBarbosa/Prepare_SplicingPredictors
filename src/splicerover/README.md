Step to generate SpliceRover finder results:

- Generate input sequences from VCF:\
`vcf2seq infile_donor.vcf.gz hg19_genome.fa out --ss donor --splicerover`\
`vcf2seq infile_acceptor.vcf.gz hg19_genome.fa out --ss acceptor --splicerover`

- Go to the [web application](http://bioit2.irc.ugent.be/rover/splicerover), update the previous generated file, select the proper model (donor or acceptor), and wait for the results. Then, save the results table in plain text and get the ready-to-annotate VCF by running:\
`python get_mutation_effect.py out_SpliceRover_donor.fa out_donor_saved_from_web.csv | bcftools sort | bgzip > splicerover_donor_to_annotate.vcf.gz`\
`python get_mutation_effect.py out_SpliceRover_acceptor.fa out_acceptor_saved_from_web.csv | bcftools sort | bgzip > splicerover_acceptor_to_annotate.vcf.gz`

- Annotate the `SpliceRover` scores into the main VCF file (to benchmark with other tools) using `vcfanno` or `bcftools annotate`
