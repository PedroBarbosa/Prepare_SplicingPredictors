Step to generate HEXplorer results:

- Generate input sequences from VCF:\
`vcf2seq infile.vcf.gz hg19_genome.fa out --hexplorer`
    - Since HEXplorer does not accept multiple fasta records, we prepare the input in a way that all variants can be tested at once. For each variant, we extract 11bp on each side, and concatenate the sequence to the previous. Since we know the exact position of each mutation (and its mapping to the original variant in VCF), it is straightforward to assign a HEXplorer score to each variant.

- The sequences (WT and mutated) should be uploaded to the [official web page](https://rna.hhu.de/HEXplorer/). Make sure that the option `Include 5nt up-/downstream` is set. **Note**: HEXplorer online interface has changed recently. That's why we have a script with `_old.py` suffix, which refers to the processing script of the output of the previous version. We have created a temporary version as `get_mutation_effects_temporary.py`, but it is likely the output format will change soon.  We will update the code as the official version gets stable.
  
- In the old setting, you would run the following command to have the proper input for the `get_mutation_effect` script:\
`cat hexplorer_results.csv | tail -n+2 | grep -v HBond > hexplorer_results_to_script.csv`
`python get_mutation_effect_old.py hexplorer_results_to_script.csv out_HEXplorer_map.tsv | bcftools sort | bgzip > HEXplorer_to_annotate.vcf.gz`\

- Now, for the temporary solution, we simply run:\
`python get_mutation_effects.py hexplorer_results.csv out_HEXplorer_map.tsv | bcftools sort | bgzip > HEXplorer_to_annotate.vcf.gz`

- Optionally, annotate the `HEXplorer` scores into the main VCF file (to benchmark with other tools) using `vcfanno` or `bcftools annotate`
