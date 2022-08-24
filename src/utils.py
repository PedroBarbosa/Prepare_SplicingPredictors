from loguru import logger
import os
import sys
from typing import Union, TextIO, Tuple
from cyvcf2 import VCF, Writer
import subprocess
import tempfile
import re
import shutil
import mygene
import pandas as pd
from pyranges import PyRanges
from pyfaidx import Fasta

PRIMARY_CHROMOSOMES = ["Y", "X", "MT", "M", "chrY", "chrX", "chrMT", "chrM"] + \
                      [str(i) for i in range(1, 23)] + \
                      ["chr" + str(i) for i in range(1, 23)]


def bed_from_vcf(vcf: Union[TextIO, str],
                 extract_genes_intervals: bool = False,
                 tx_field_name: str = "Feature",
                 gene_field_name: str = "Gene",
                 genome: str = "hg38"):
    """
    Extracts bed intervals where variants occur. By default,
    output represents exact intervals spanned by the variant
    and strand information is retrieved from VEP annotations,
    if they exist (STRAND attribute within `ANN` or `CSQ` INFO
    subfield) using the first VEP consequence.


    If `extract_genes_intervals` is `True`, full gene-length
    intervals will be retrieved. gene IDs will be obtained from VEP
    annotations (if they exist) from first consequence described.
    From there, gene information will be fetched using mygene package
    so that strand aware annotations are properly produced. If
    geneID is not found for a given variant (e.g. intergenic),
    variant is skipped.

    :param Union[TextIO, str] vcf: Path to the VCF file
    :param bool extract_genes_intervals: Whether full gene-length
        intervals where variant occurs should be retrieved.
        Default: `False`, exact intervals are returned.
    :param str tx_field_name: Name of the field that
        corresponds to the transcript ID within VEP
        annotations. Default: `Feature`
    :param str gene_field_name: Name of the field that
        corresponds to the gene ID within VEP annotations
        and used when `extract_genes_intervals`
        is `True`.Default: `Gene`
    :param str genome: Human genome build to fetch
        genomic coordinates. Default: "hg38" if
        `extract_genes_intervals` is `True`.

    :return Tuple: 2 elements tuple with pyranges
        object with bed (0-based) intervals (1st
        element) and a list containing the unique
        transcript IDs affected by variants, considering
        the first consequence of the VEP annotations, if
        they exist
    """

    logger.info("Extracting variant intervals")
    assert genome in ['hg38', 'hg19'], "Invalid genome build set."
    vep_markers = ['ANN', 'CSQ']
    vcf_data = VCF(vcf)

    if any(vcf_data.contains(x) for x in vep_markers):
        for field in vcf_data.header_iter():
            if field["HeaderType"] == "INFO" and field["ID"] in vep_markers:
                vep_annotations = field["Description"].split("Format:")[1][:-1].strip().split("|")

                assert tx_field_name in vep_annotations, "Make sure 'tx_field_name' is present in " \
                                                         "VEP annotations. '{}' was set, but it " \
                                                         "doesn't exist in the VCF header.".format(tx_field_name)
                tx_id_idx = vep_annotations.index(tx_field_name)

                # Extract full gene length intervals
                if extract_genes_intervals:
                    assert gene_field_name in vep_annotations, "Make sure 'gene_field_name' is present in " \
                                                               "VEP annotations. '{}' was set, but it doesn't " \
                                                               "exist in the VCF header.".format(gene_field_name)
                    gene_id_idx = vep_annotations.index(gene_field_name)
                    return _get_full_gene_length_intervals(vcf_data, tx_id_idx, gene_id_idx, genome)

                # Extract exact variant intervals
                else:
                    assert "STRAND" in vep_annotations, "'STRAND' attribute must exist in VEP annotations"
                    try:
                        indexes = (vep_annotations.index("STRAND"),
                                   vep_annotations.index("Gene"),
                                   vep_annotations.index("SYMBOL"),
                                   vep_annotations.index("Feature"),
                                   vep_annotations.index("HGVSc"),
                                   vep_annotations.index("EXON"),
                                   vep_annotations.index("INTRON"),
                                   vep_annotations.index("VARIANT_CLASS"))

                    # If one of the fields doesn't exist, just use STRAND
                    except ValueError:
                        indexes = vep_annotations.index("STRAND")

                    return _get_variant_intervals(vcf_data, tx_id_idx=tx_id_idx, target_fields_idx=indexes)

    else:
        logger.info("VEP annotations were not found, strand "
                     "information won't be retrieved. Forward "
                     "strand will be extracted for all the "
                     "variants and `extract_genes_intervals`"
                     "is disabled.")
        return _get_variant_intervals(vcf_data)


def _get_variant_intervals(vcf_data: VCF,
                           tx_id_idx: int = None,
                           target_fields_idx: Union[Tuple, int] = None):
    """
    Given a cyvcf2.VCF object, returns
    bed intervals (works for SNV and Indels)
    :param VCF vcf_data: cyvcf2.VCF object containing
        the variants to iterate over.
    :param int tx_id_idx: Index where transcriptIDs
        are present within the VEP annotation field.
        Default: None, assumes no VEP annotations
        are available
    :param Tuple target_fields_idx: Indexes where strand
        info (and additional fields) are within the VEP
        annotation fields. Default: None, assumes no
        VEP annotations are available
    :return PyRanges: pyranges object of genomic intervals
        representing input variants
    :return list: List containing the unique transcript
        IDs affected by variants, considering the
        first consequence of the VEP annotations, if
        they exist
    """
    transcript_ids = set()
    strand_map = {"-1": "-", "1": "+"}
    strand_idx = None
    intervals = []

    #TODO Need to check variants with no genic consequences within VEP
    # Perhaps create a way to check further consequences or just remove
    # those variants. (e.g. gnomAD benign variants in target genes)

    # If no VEP annotation, strand won't be retrieved
    cols = ['Chromosome', 'Start', 'End', 'Pos', 'ID', 'Score', 'Strand']

    # If strand and additional target fields exist
    if isinstance(target_fields_idx, Tuple):
        strand_idx = target_fields_idx[0]
        cols = ['Chromosome', 'Start', 'End', 'Pos', 'ID',
                'HGVSc', 'Strand', 'gene_id', 'gene_name',
                'transcript_id', 'Exon_number',
                'Intron_number', 'Type']

    cols.extend(['REF', 'ALT'])
    for record in vcf_data:

        ann = record.INFO.get("ANN") if vcf_data.contains("ANN") else record.INFO.get("CSQ")
        if ann:
            ann_fields = ann.split(",")[0].split("|")
            tx_id = ann_fields[tx_id_idx]
            transcript_ids.add(tx_id)
            if target_fields_idx is not None:
                strand = ann_fields[strand_idx] if ann else None
        else:
            strand = None

        if len(record.ALT) > 1:
            raise ValueError('Multiallelic variant ({}, {}) found in '
                             'the VCF file. Split multiallelic variants '
                             'before extracting bed intervals from a VCF file.')
        record_out = []
        ALT = record.ALT[0]

        # Deletions: Output interval comprises the
        # length of the reference allele
        if len(record.REF) > len(ALT):
            record_out.extend([record.CHROM, record.POS - 1, record.POS + len(record.REF) - 1,
                               record.POS, record.ID])

        # SNVs, MNPs and Insertions: Output interval
        # comprises the length of additional bases
        # (if insertion) in the alt allele
        else:
            record_out.extend([record.CHROM, record.POS - 1, record.POS + len(record.ALT[0]) - 1,
                               record.POS, record.ID])

        if strand:
            # If additional fields exist in VEP annotations
            if isinstance(target_fields_idx, Tuple):
                gene_id = ann_fields[target_fields_idx[1]]
                gene_name = ann_fields[target_fields_idx[2]]
                tx_id = ann_fields[target_fields_idx[3]]
                hgvs = ann_fields[target_fields_idx[4]]
                exon = ann_fields[target_fields_idx[5]]
                intron = ann_fields[target_fields_idx[6]]
                type = ann_fields[target_fields_idx[7]]
                record_out.extend([hgvs, strand_map[strand], gene_id, gene_name, tx_id, exon, intron, type])
            else:
                record_out.extend(["0", strand_map[strand]])
        else:
            record_out.extend(["0", "+"])
        record_out.extend([record.REF, ALT])
        intervals.append(record_out)

    transcript_ids = None if len(transcript_ids) == 0 else list(transcript_ids)
    return PyRanges(pd.DataFrame.from_records(intervals, columns=cols)).sort(), transcript_ids


def annotate_vcf(vcf: str,
                 df: Union[pd.DataFrame, PyRanges],
                 col_to_add: str,
                 annot_name: str = None,
                 description: str = None,
                 _type: str = 'String',
                 _number: str = '1',
                 output_file: str = "annotated.vcf"
                 ):
    """
    Annotates a VCF with a field present
    in the provided df and writes new file.
    Df must contain REF and ALT alleles so
    that the same exact variants are compared.

    :param str vcf: Path to the VCF file
    :param Union[pd.DataFrame, PyRanges] df: Input
    df with field to add.
    :param str col_to_add: Column from `df` to
    add to VCF
    :param str annot_name: Name of new annotation
    in VCF. Default: `None`, name will be the same
    as `col_to_add`.
    :param str description: Description of the new
    annotation to add to the VCF header. Default: `.`
    :param str _type: Type for the new INFO field.
        Default: 'String'. Possible values:
        [Integer, Float, Flag, Character, String]
    :param str _number: Number of values that can be
        included with the INFO field. Default: `1`,
        values in `col_to_add` refer to a single element.
    :param str output_file: Name of the output file.
    Default: `annotated.vcf
    :return:
    """
    # Need to update method to annotate
    # multiple fields at once, if given
    required_cols = ['Chromosome', 'Pos', 'REF', 'ALT']
    assert all(x in df.columns for x in required_cols), 'Df does not contain ' \
                                                        'required cols: {}\n'.format(required_cols)

    assert col_to_add in df.columns, "Column {} must exist in input df.".format(col_to_add)
    assert _type in ["Integer", "Float", "Flag", "Character", "String"], "'_type' must be a valid value " \
                                                                         "from VCF specifications."

    new_annotation = annot_name if annot_name else col_to_add

    vcf_data = VCF(vcf)

    # Ensure that new annotation is not yet in VCF
    existing_annotations = []
    for field in vcf_data.header_iter():

        # if VEP annotations
        if field["HeaderType"] == "INFO" and (field["ID"] == "ANN" or field["ID"] == "CSQ"):
            vep_subfields = field["Description"].split("Format:")[1][:-1].strip().split("|")
            existing_annotations.extend(vep_subfields)

        elif field["HeaderType"] == "INFO":
            existing_annotations.append(field['ID'])

    assert new_annotation not in existing_annotations, "New annotation name ({}) already " \
                                                       "exists in VCF. Please set a different " \
                                                       "name in 'col_to_add' or 'annot_name' " \
                                                       "argument. List of existing fields in" \
                                                       "VCF:\n{}".format(new_annotation, existing_annotations)

    # Add new info to header
    vcf_data.add_info_to_header({'ID': new_annotation,
                                 'Description': description,
                                 'Type': _type,
                                 'Number': _number})

    # Setup output writer
    w = Writer(output_file, vcf_data)

    # Iterate over VCF
    for record in vcf_data:

        v_id = record.CHROM + "_" + str(record.POS) + "_" + record.REF + "_" + record.ALT[0]
        query = df[(df.Chromosome == record.CHROM) &
                   (df.Pos == record.POS) &
                   (df.REF == record.REF) &
                   (df.ALT == record.ALT[0])]

        if query.empty:
            record.INFO[new_annotation] = "."

        else:
            _value = query[col_to_add]

            # multiple queries found
            if _value.size > 1:
                raise ValueError('More than one entry in df exists to annotate '
                                 'the variant ({}). Not sure which one to use. '
                                 'Values at the target column: {}'.format(v_id, _value.tolist()))

            elif _value.size == 1:
                record.INFO[new_annotation] = _value.iloc[0]

            else:
                raise ValueError("Unknown problem found for variant id {} "
                                 "and queried value in df ({})".format(v_id, query))

        w.write_record(record)
    vcf_data.close()
    w.close()


def _get_full_gene_length_intervals(vcf_data: VCF, tx_id_idx: int, gene_id_idx: int, genome: str):
    """
    Returns unique full gene-length intervals
    spanning the variants present in a `cyvcf2.VCF`
    object

    :param VCF vcf_data: cyvcf2.VCF object containing
        the variants to iterate over.
    :param int tx_id_idx: Index where transcriptIDs
        are present within the VEP annotation field
    :param int gene_id_idx: Index where geneIDs are
        present within the VEP annotation field
    :param str genome: Genome build to fetch data
    :return PyRanges: pyranges object with gene
        coordinates that overlap the input variants
    :return list: List containing the unique transcript
        IDs affected by variants, considering the
        first consequence of the VEP annotations
    """
    gene_ids, transcript_ids = set(), set()
    for record in vcf_data:
        annotation = record.INFO.get("ANN") if vcf_data.contains("ANN") else record.INFO.get("CSQ")
        if annotation:
            g_id = annotation.split(",")[0].split("|")[gene_id_idx]
            tx_id = annotation.split(",")[0].split("|")[tx_id_idx]
            gene_ids.add(g_id)
            transcript_ids.add(tx_id)

    # Generate bed entries from geneIDs
    mg = mygene.MyGeneInfo()
    fields = ["genomic_pos.chr", "genomic_pos.start", "genomic_pos.end", "symbol", "genomic_pos.strand"]
    if genome == "hg19":
        fields = [re.sub(r'_pos', '_pos_hg19', x) for x in fields]
    ens_map = mg.getgenes(ids=gene_ids,
                          fields=fields,
                          returnall=True,
                          as_dataframe=True,
                          size=1,
                          species="human")
    cols = ["Chromosome", "Start", "End", "Name", "Score", "Strand"]
    to_bed = ens_map.reset_index().rename(columns={fields[0]: 'Chromosome',
                                                   fields[1]: 'Start',
                                                   fields[2]: 'End',
                                                   fields[3]: 'Name',
                                                   fields[4]: 'Strand',
                                                   'query': "Score"})[cols]

    # Tricky IDs
    tricky_ids = to_bed[to_bed.isnull().any(axis=1)]['Score']
    d = ens_map.loc[tricky_ids][['genomic_pos', 'symbol']].to_dict()

    _to_df = {}
    for gene_id, matched in d['genomic_pos'].items():
        for m in matched:
            if m['chr'] in PRIMARY_CHROMOSOMES:
                m['symbol'] = d['symbol'][gene_id]
                _to_df[gene_id] = m

    _df = pd.DataFrame.from_dict(_to_df, orient='index').reset_index().rename(
        columns={'chr': 'Chromosome',
                 'start': 'Start',
                 'end': 'End',
                 'symbol': 'Name',
                 'index': 'Score',
                 'strand': 'Strand'})

    # Merge
    to_bed = to_bed[~to_bed.isnull().any(axis=1)]
    to_bed[["Start", "End", "Strand"]] = to_bed[["Start", "End", "Strand"]].astype(int)
    to_bed = pd.concat([to_bed, _df])
    to_bed['Chromosome'] = 'chr' + to_bed['Chromosome']
    to_bed['Start'] -= 1
    to_bed['Strand'] = to_bed['Strand'].replace({-1: "-", 1: "+"})
    return PyRanges(to_bed).sort(), list(transcript_ids)



def _fix_chr(vcf_object: VCF, out_dir: str, add: bool = True):
    """
    Fix chromosome naming in a VCF file

    :param VCF vcf_object: cyvcf2 VCF object
    :param str out_dir: Directory to write the new VCF,
        in case a fix on chromosome names was required.
    :param bool add: Whether 'chr' string should be added.
        Default: True
    :return:
    """
    all_chrom = vcf_object.seqnames
    header = vcf_object.raw_header
    new_header = []
    contig_id_found = False
    for row in header.split("\n"):
        if not contig_id_found and row.startswith("##contig=<ID="):
            contig_id_found = True

        if row.startswith("##contig=<ID=") and row.split("##contig=<ID=")[1][:-1] in PRIMARY_CHROMOSOMES:

            chrom = row.split("##contig=<ID=")[1][:-1]
            if chrom.endswith("MT"):
                logger.info("MT chr found in VCF file. It will be replaced with M.")
                chrom = chrom[:-1]

            if chrom in all_chrom or "chr" + chrom in all_chrom:
                new_header.append("##contig=<ID=chr" + chrom + ">") if add else \
                    new_header.append("##contig=<ID=chr" + re.sub(r'^chr', '', chrom) + ">")

        # If contig names not in header at all
        elif row.startswith("#CHROM") and not contig_id_found:
            for c in PRIMARY_CHROMOSOMES:
                if add and c.startswith("chr") and c != "chrMT":
                    new_header.append("##contig=<ID=" + c + ">")
                elif not add and not c.startswith("chr") and c != "MT":
                    new_header.append("##contig=<ID=" + c + ">")
            new_header.append(row)

        else:
            new_header.append(row)

    with tempfile.NamedTemporaryFile() as tmp:
        tmp.write('\n'.join(new_header).encode("UTF-8"))
        for record in vcf_object:
            line = str(record).split()
            if line[0] in PRIMARY_CHROMOSOMES:
                line[0] = line[0][:-1] if line[0].endswith("MT") else line[0]
                line[0] = "chr" + line[0] if add else re.sub(r'^chr', '', line[0])

            tmp.write(('\t'.join(line) + '\n').encode("UTF-8"))

        fn = os.path.join(out_dir, "0_variants_chr_fixed.vcf")
        shutil.copy(tmp.name, fn)
    tmp.close()
    subprocess.run(["bgzip", "-f", fn])
    subprocess.run(["tabix", "-pvcf", fn + ".gz"])
    return fn + ".gz"
