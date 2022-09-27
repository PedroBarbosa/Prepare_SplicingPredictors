Step to generate ESEfinder results:

- Generate input sequences from VCF (if ):\
`vcf2seq infile.vcf.gz hg19_genome.fa out --esefinder`
    - If there are many variants in VCF, it's better to split them in several files, since the ESEfinder server has a limit in the maximum number of nucleotides (all sequences combined) allowed as input.
    - The input is prepared in a way that all variants can be tested at once. For each variant, we extract 10bp on each side, and ESEfinder will score at each position of the input the value for SR protein.

- The sequences should be uploaded to the [official web page](http://krainer01.cshl.edu/cgi-bin/tools/ESE3/esefinder.cgi). Remove any empty line from the FASTA file (e.g. last line) and select two options in the `Output information` section: `Report all scores in each sequence (instead of hits above the thresholds)` and `Output as a plain text file`. Copy the output from clipboard to a new file and run the `get_mutation_effect` script:\
`python get_mutation_effect.py ESEfinder_out.tsv`

- Optionally, annotate the `ESEfinder` scores into the main VCF file (to benchmark with other tools) using `vcfanno` or `bcftools annotate`
