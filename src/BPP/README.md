Step to generate BPP results:

- Generate input sequences from VCF (500bp upstream of acceptor site):\
`vcf2seq infile.vcf.gz hg19_genome.fa out --bpp`

- Download the [official repo](https://github.com/zhqingit/BPP), and run the model accordingly:\
`BP_PPT.py -b demo/pwmBP_human.txt -p demo/scPPT_human.txt -i out_BPP.fa > BPP_out.tsv`

- Process output and automatically generate a VCF :\
`python get_mutation_effect.py BPP_out.tsv | bcftools sort | bgzip > BPP.vcf.gz`

- Optionally, annotate the `BPP` scores into the main VCF file (to benchmark with other tools) using `vcfanno` or `bcftools annotate`
