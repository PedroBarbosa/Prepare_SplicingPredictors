import os 
from cyvcf2 import VCF
from pyfaidx import Fasta

PRIMARY_CHROMOSOMES = ["Y", "X", "MT", "M", "chrY", "chrX", "chrMT", "chrM"] + \
                      [str(i) for i in range(1, 23)] + \
                      ["chr" + str(i) for i in range(1, 23)]
                      
def validate_vep(vcf_f: str,
                 fasta: Fasta = None):
    """
    Reads a VCF, validates it and fixes chromosome names if necessary.
    By default, with `Fasta` = `None`, 'chr' string will be excluded
    from the VCF, in case it exist.

    :param str vcf_f: VCF file

    :return list: Returns indexes of VEP annotations required 
    """

    if not os.path.isfile(vcf_f):
        raise FileNotFoundError("VCF file does not exist")

    vep_markers = ['ANN', 'CSQ']
    vcf_data = VCF(vcf_f)

    if any(vcf_data.contains(x) for x in vep_markers):
        for field in vcf_data.header_iter():
            if field["HeaderType"] == "INFO" and field["ID"] in vep_markers:
                vep_annotations = field["Description"].split("Format:")[1][:-1].strip().split("|")

                for attr in ['STRAND', 'HGVSc', 'EXON', 'INTRON', 'Feature']:
                    assert attr in vep_annotations, "{} attribute must exist in VEP annotations".format(attr)
                
                indexes = {'strand': vep_annotations.index("STRAND"),
                           'hgvsc': vep_annotations.index("HGVSc"),
                           'exon': vep_annotations.index("EXON"),
                           'intron': vep_annotations.index("INTRON"),
                           'tx_id': vep_annotations.index("Feature")}

    else:
        raise ValueError("VEP annotations were not found. Please run ensembl VEP before. Don't forget to set the '--hgvs' flag")

    return indexes