Step to generate HEXplorer results:

- Generate input sequences from VCF:\
`vcf2seq infile.vcf.gz hg19_genome.fa out --hexplorer`
    - Since HEXplorer does not accept multiple fasta records, we prepare the input in a way that all variants can be tested at once. For each variant, we extract 60bp on each side so that the length of each input sequence sums to the median length of a coding human exon length (120bp), and concatenate the sequence to the previous. Since we know the exact position of each mutation (and its mapping to the original variant in VCF), it is straightforward to assign a HEXplorer score to each variant.

- The sequences (WT and mutated) should be uploaded to the [official web page](https://www2.hhu.de/rna/html/hexplorer_score.php) (**Note**: HEXplorer online interface has changed recently. We did not test the new version yet). Calculate the scores and download the results.
  
  
-  Then, run the following command to have the proper input for the `get_mutation_effect` script:\
`cat hexplorer_results.csv | tail -n+2 | grep -v HBond > hexplorer_results_to_script.csv`
`python get_mutation_effect hexplorer_results_to_script.csv out_HEXplorer_map.tsv | bcftools sort | bgzip > HEXplorer_to_annotate.vcf.gz`\

- Optionally, annotate the `HEXplorer` scores into the main VCF file (to benchmark with other tools) using `vcfanno` or `bcftools annotate`
