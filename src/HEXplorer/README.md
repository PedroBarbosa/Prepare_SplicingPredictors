Step to generate HEXplorer results:

- Generate input sequences from VCF:\
`vcf2seq infile.vcf.gz hg19_genome.fa test --hexplorer`
    - Since HEXplorer does not accept multiple fasta records, we prepare the input in a way that multiple variants can be tested at once. For each variant, we extract 11bp on each side, and concatenate the sequence to the previous. Since we know the exact position of each mutation (and its mapping to the original variant in VCF), it is straightforward to assign a HEXplorer score to each variant.

- The sequences (WT and mutated) should be uploaded to the [official web page](https://rna.hhu.de/HEXplorer/). Make sure that the option `Include 5nt up-/downstream` is set.

    **Note1**: HEXplorer online interface has changed recently. That's why we have a script with `_old.py` and `_temporary.py` suffixes, which refers to the processing script of the output from previous versions.

    **Note2**: Because the server has limits regarding the input sequence length, we split the input sequences if there are more than 50 variants to analyze.

- After downloading the results table, simply run:\
`python get_mutation_effects.py hexplorer_results.csv test_HEXplorer_map.tsv | bcftools sort | bgzip > HEXplorer_to_annotate.vcf.gz`. 

    This script will assign a prediction value to each variant by summing the HZEI differences betwen the mutated - WT of all hexamers overlapping the variant.

- Finally, annotate the `HEXplorer` scores into the main VCF file (to benchmark with other tools) using `vcfanno` or `bcftools annotate`

------------

- In the old setting, we would return the max HZEI difference (instead of the sum) found in an hexamer overlapping the variant: 
`cat hexplorer_results.csv | tail -n+2 | grep -v HBond > hexplorer_results_to_script.csv`
`python get_mutation_effect_old.py hexplorer_results_to_script.csv test_HEXplorer_map.tsv | bcftools sort | bgzip > HEXplorer_to_annotate.vcf.gz`\

- Similarly, for the temporary solution, we would run:\
`python get_mutation_effects_temporary.py hexplorer_results.csv test_HEXplorer_map.tsv | bcftools sort | bgzip > HEXplorer_to_annotate.vcf.gz`