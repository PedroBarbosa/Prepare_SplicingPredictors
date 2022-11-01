# Preparing splicing-related models
A simple utility to process VCF files and generate correct input for a set of models that require sequence information (or that require specific format to run a web-based application).

### Motivation
There are a lot of available methods that predict splicing-related information (e.g. splice sites, branchpoints, splicing regulatory elements). Since they were not originally designed to predict the effect of genetic variants, it is not straightforward to use these models for that task. This tool simplifies that goal: generates reference and mutated sequences from VCF files in the proper format to run several models for all variants at once, and contains utilities to process the output and generate a VCF with a final score (usually mutated allele - reference allele).

### Requirements 
The variants should be annotated with ensembl VEP so that strand information can be retrieved (and therefore the proper sequence context of the variant can be extracted).

### Instalation
```
git clone https://github.com/PedroBarbosa/Prepare_SplicingPredictors.git
cd Prepare_SplicingPredictors
conda env create --file conda_environment.yaml 
conda activate prepareSplicingTools
pip install .
```
### Running 
To run this utility, just call the `vcf2seq` and select the models you want to generate input for (check the available options with `vcf2seq --help`).
```
vcf2seq input.vcf.gz reference_genome.fa outbasename --maxentscan --splicerover ...
```

For models that predict splice sites, it may be necessary to set the splice site flag (`--ss donor`, `--ss acceptor`).
Then, within each model folder (`src` folder in this repo), there are instructions on how to run each model and a script (`get_mutation_effects.py`) to process the output and generate a VCF with the predictions.

Note: Do not change the fasta headers of the generated sequences, since the `get_mutation_effects.py` scripts require original names for proper processing.

### Supported models
#### General methods
* [regSNP-intron](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-019-1847-4)
  
#### Splice site prediction
* [Splice2deep](https://www.sciencedirect.com/science/article/pii/S2590158320300097)
* [SpliceRover](https://pubmed.ncbi.nlm.nih.gov/29931149/)
* [DSSP](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7265986/)
* [Spliceator](https://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-021-04471-3)
* [MaxEntScan](https://github.com/kepbod/maxentpy)

#### Splicing regulatory elements
* [ESEfinder](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC169022/)
* [ESRseq](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3149502/)
* [HEXplorer](https://pubmed.ncbi.nlm.nih.gov/25147205/)

#### Branchpoint signals
* [SVM-BPfinder](https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1001016)
* [BPP](https://academic.oup.com/bioinformatics/article/33/20/3166/3870482)
* [BPHunter](https://www.pnas.org/doi/10.1073/pnas.2211194119)
* 
### Limitations
Only single-nucleotide variants (SNVs) are supported.
