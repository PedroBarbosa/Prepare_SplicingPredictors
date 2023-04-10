Software can be run using a [web interface](https://hgidsoft.rockefeller.edu/BPHunter/index.php), but we show here how to run it locally:

- Clone [git repository](https://github.com/casanova-lab/BPHunter.git) and make sure `bedtools` is installed.

- Download reference datasets from [here](https://hgidsoft.rockefeller.edu/BPHunter/standalone.html) and put them in the repository folder.

- Unzip VCF file (e.g. with `zcat`)

- `python BPHunter_VCF.py -i input.vcf -g GRCh37 -t all`

- Generate a VCF from the raw BPHunter output:\
`python get_mutation_effect.py input_BPHunter_output_1681231.txt | bcftools sort | bgzip > BPHunter.vcf.gz`

Alternatively, if you want to score all variants in the initial VCF file (and thus evaluate its performance fairly), you can shift variants with 0 score to 1 (by providing `shift` as the 2nd argument of the script), and then annotate the initial VCF with `vcfanno` by including a postannotation step to annotate all unnanotated variants with 0. This shifting procedure was personally suggested by BPHunter author:

`python get_mutation_effect.py input_BPHunter_output_1681231.txt shift | bcftools sort | bgzip > BPHunter.vcf.gz`

`vcfanno -lua function.lua vcfanno.conf input.vcf.gz | bgzip > BPHunter_fully_annotated.vcf`
