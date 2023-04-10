import sys

from vcf_utils import *
import cyvcf2
import hgvs.parser
from Bio.Seq import Seq
from pyfaidx import Fasta
import argparse
import logging

logging.basicConfig(
    stream=sys.stdout, level=logging.INFO, format="%(asctime)s %(message)s"
)


def _no_ss(id, seq_wt, seq_mut, ss, tool):
    """
    Models that test the value at each position (such as in SpliceAI),
    don't require specific ss motif at the 4-mer (e.g. Spliceator).
    """
    if ss == "donor":
        _ss = "GT"
    else:
        _ss = "AG"

    if all(_ss not in x for x in [seq_wt, seq_mut]):
        if tool in ["Splice2Deep", "DSSP"]:
            logging.info(
                "[{} {}] Variant {} discarded. Reason: no {} dinucleotide formed at the 4-mer around the mutation (wt={}, mut={}).".format(
                    tool, ss, id, _ss, seq_wt, seq_mut
                )
            )
            return True
        elif tool in ["Spliceator", "SpliceRover"]:
            logging.info(
                "[{} {}] Warning variant {}: No {} dinucleotide formed at the 4-mer around the mutation (wt={}, mut={}).".format(
                    tool, ss, id, _ss, seq_wt, seq_mut
                )
            )
            return False
    else:
        return False


def _dict_to_fasta(d: dict, outfile: str):
    """
    Writes a dict (keys are the header, values are the seqs)
    to a fasta file
    :param dict d: Input dictionary with sequences
    :param str outfile: Output file
    """
    if d and len(d) > 0:
        f = open(outfile, "w")
        for k, v in d.items():
            f.write(">{}\n{}\n".format(k, v))


class DataProcessing(object):
    def __init__(
        self,
        vcf: str,
        vep_indexes: dict,
        vep_tag: str,
        fasta: Fasta,
        outbasename: str,
        ss: str = None,
        svm_bp_finder: bool = False,
        bpp: bool = False,
        maxentscan: bool = False,
        spliceator: bool = False,
        splice2deep: bool = False,
        splicerover: bool = False,
        dssp: bool = False,
        hexplorer: bool = False,
        esefinder: bool = False,
        esrseq: bool = False,
        surrounding: int = None,
    ):
        self.hp = hgvs.parser.Parser()
        self.vep_tag = vep_tag
        self.vep_indexes = vep_indexes
        self.fasta = Fasta(fasta)
        self.outbasename = outbasename
        self.ss = ss
        self.svm_bp_finder = svm_bp_finder
        self.bpp = bpp
        self.maxentscan = maxentscan
        self.spliceator = spliceator
        self.splice2deep = splice2deep
        self.dssp = dssp
        self.splicerover = splicerover
        self.hexplorer = hexplorer
        self.esefinder = esefinder
        self.esrseq = esrseq
        self.surrounding = surrounding
        self._iterate_vcf(vcf)

    def _iterate_vcf(self, vcf):
        (
            fasta_bpp,
            fasta_svm_bp,
            fasta_maxentscan,
            fasta_spliceator,
            fasta_splice2deep,
            fasta_dssp,
            fasta_splicerover,
            fasta_esefinder,
            fasta_esrseq,
            fasta_surrounding,
        ) = ({}, {}, {}, {}, {}, {}, {}, {}, {}, {})
        hexplorer_mapping, hexplorer_ref, hexplorer_mut = {}, "", ""

        for record in VCF(vcf):
            if len(record.ALT) > 1:
                raise ValueError(
                    "Multi allelic records are not allowed. Please split them before."
                )
            vep_annotation = record.INFO.get(self.vep_tag)

            if vep_annotation:
                # First consequence is selected
                vep_annotation = vep_annotation.split(",")[0]
            else:
                continue

            if self.svm_bp_finder:
                fasta_svm_bp = self.BranchpointTools(
                    record, vep_annotation, out_dict=fasta_svm_bp, tool="SVM_BP_finder"
                )

            if self.bpp:
                fasta_bpp = self.BranchpointTools(
                    record, vep_annotation, out_dict=fasta_bpp, tool="BPP"
                )

            if self.maxentscan:
                fasta_maxentscan = self.MaxEntScan(
                    record, vep_annotation, out_dict=fasta_maxentscan
                )
            if self.spliceator:
                fasta_spliceator = self.Spliceator(
                    record, vep_annotation, out_dict=fasta_spliceator
                )

            if self.splice2deep:
                fasta_splice2deep = self.Splice2Deep(
                    record, vep_annotation, out_dict=fasta_splice2deep
                )

            if self.dssp:
                fasta_dssp = self.DSSP(record, vep_annotation, out_dict=fasta_dssp)

            if self.splicerover:
                fasta_splicerover = self.SpliceRover(
                    record, vep_annotation, out_dict=fasta_splicerover
                )

            if self.esefinder:
                fasta_splicerover = self.ESEfinder(
                    record, vep_annotation, out_dict=fasta_esefinder
                )

            if self.hexplorer:
                hexplorer_mapping, hexplorer_ref, hexplorer_mut = self.HEXplorer(
                    record,
                    vep_annotation,
                    variant_seq_map=hexplorer_mapping,
                    hexplorer_ref_seq=hexplorer_ref,
                    hexplorer_mut_seq=hexplorer_mut,
                )

            if self.esrseq:
                fasta_esrseq = self.ESRseq(
                    record, vep_annotation, out_dict=fasta_esrseq
                )

            if self.surrounding is not None:
                assert self.surrounding > 0, "--surrounding must be a positive integer"
                fasta_surrounding = self.extract_surrounding(
                    record, vep_annotation, out_dict=fasta_surrounding
                )

        ###################
        ## WRITE OUTPUT  ##
        ###################
        if self.svm_bp_finder:
            _dict_to_fasta(fasta_svm_bp, "{}_SVM_BP_finder.fa".format(self.outbasename))

        if self.bpp:
            _dict_to_fasta(fasta_bpp, "{}_BPP.fa".format(self.outbasename))

        if self.maxentscan:
            _dict_to_fasta(
                fasta_maxentscan,
                "{}_MaxEntScan_{}.fa".format(self.outbasename, self.ss),
            )

        if self.spliceator:
            _dict_to_fasta(
                fasta_spliceator,
                "{}_Spliceator_{}.fa".format(self.outbasename, self.ss),
            )

        if self.splice2deep:
            _dict_to_fasta(
                fasta_splice2deep,
                "{}_Splice2Deep_{}.fa".format(self.outbasename, self.ss),
            )

        if self.dssp:
            _dict_to_fasta(
                fasta_dssp, "{}_DSSP_{}.fa".format(self.outbasename, self.ss)
            )

        if self.splicerover:
            _dict_to_fasta(
                fasta_splicerover,
                "{}_SpliceRover_{}.fa".format(self.outbasename, self.ss),
            )

        if self.esefinder:
            _dict_to_fasta(fasta_esefinder, "{}_ESEfinder.fa".format(self.outbasename))

        if self.esrseq:
            _dict_to_fasta(fasta_esrseq, "{}_ESRseq.fa".format(self.outbasename))

        if self.hexplorer:
            f = open("{}_HEXplorer_map.tsv".format(self.outbasename), "w")
            for v_id, seq_mut_pos in hexplorer_mapping.items():
                f.write("{}\t{}\n".format(v_id, seq_mut_pos))
            f.close()
            f = open("{}_HEXplorer_input.txt".format(self.outbasename), "w")
            f.write(
                ">seqs_{}\n{}\n{}".format(
                    self.outbasename, hexplorer_ref, hexplorer_mut
                )
            )
            f.close()

        if self.surrounding is not None:
            _dict_to_fasta(
                fasta_surrounding, "{}_surrounding_{}bp.fa".format(self.outbasename, self.surrounding)
            )
            
    def extract_surrounding(
        self, record: cyvcf2.Variant, vep_annotation: str, out_dict: dict
    ):
        """
        Sequences of surrounding * 2 + 1 will be generated,
        with the variant located in the middle position.
        """

        id = (
            record.CHROM
            + "_"
            + str(record.POS)
            + "_"
            + record.REF
            + "_"
            + record.ALT[0]
        )
        header_wt = id + "_WT"
        header_mut = id + "_Mutated"
        strand = vep_annotation.split("|")[self.vep_indexes["strand"]]
        if record.var_type in ["snp"]:
            if strand == "1":
                seq_wt = str(
                    self.fasta[record.CHROM][
                        record.POS
                        - self.surrounding
                        - 1 : record.POS
                        + self.surrounding
                    ]
                )
                seq_mut = (
                    seq_wt[: self.surrounding]
                    + record.ALT[0]
                    + seq_wt[self.surrounding + 1 :]
                )

            elif strand == "-1":
                seq_wt = str(
                    -self.fasta[record.CHROM][
                        record.POS
                        - self.surrounding
                        - 1 : record.POS
                        + self.surrounding
                    ]
                )
                seq_mut = (
                    seq_wt[: self.surrounding]
                    + str(Seq(record.ALT[0]).complement())
                    + seq_wt[self.surrounding + 1 :]
                )

            out_dict[header_wt] = seq_wt
            out_dict[header_mut] = seq_mut

        else:
            logging.info(
                "[Surrounding] Variant {} discarded. Reason: There is only support for SNVs".format(
                    id
                )
            )
        return out_dict

    def ESRseq(self, record: cyvcf2.Variant, vep_annotation: str, out_dict: dict):
        """
        11bp sequences will be generated with each
        variant being located in the middle position (6).
        """
        id = (
            record.CHROM
            + "_"
            + str(record.POS)
            + "_"
            + record.REF
            + "_"
            + record.ALT[0]
        )
        header_wt = id + "_WT"
        header_mut = id + "_Mutated"
        strand = vep_annotation.split("|")[self.vep_indexes["strand"]]
        if record.var_type in ["snp"]:
            if strand == "1":
                seq_wt = str(self.fasta[record.CHROM][record.POS - 6 : record.POS + 5])
                seq_mut = seq_wt[:5] + record.ALT[0] + seq_wt[6:]

            elif strand == "-1":
                seq_wt = str(-self.fasta[record.CHROM][record.POS - 6 : record.POS + 5])
                seq_mut = seq_wt[:5] + str(Seq(record.ALT[0]).complement()) + seq_wt[6:]

            out_dict[header_wt] = seq_wt
            out_dict[header_mut] = seq_mut

        else:
            logging.info(
                "[ESRseq] Variant {} discarded. Reason: There is only support for SNVs".format(
                    id
                )
            )
        return out_dict

    def ESEfinder(self, record: cyvcf2.Variant, vep_annotation: str, out_dict: dict):
        """
        21bp sequences will be generated with each
        variant being located in the middle position (61).
        """
        id = (
            record.CHROM
            + "_"
            + str(record.POS)
            + "_"
            + record.REF
            + "_"
            + record.ALT[0]
        )
        header_wt = id + "_WT"
        header_mut = id + "_Mutated"
        strand = vep_annotation.split("|")[self.vep_indexes["strand"]]
        if record.var_type in ["snp"]:
            if strand == "1":
                seq_wt = str(
                    self.fasta[record.CHROM][record.POS - 11 : record.POS + 10]
                )
                seq_mut = seq_wt[:10] + record.ALT[0] + seq_wt[11:]

            elif strand == "-1":
                seq_wt = str(
                    -self.fasta[record.CHROM][record.POS - 11 : record.POS + 10]
                )
                seq_mut = (
                    seq_wt[:10] + str(Seq(record.ALT[0]).complement()) + seq_wt[11:]
                )

            out_dict[header_wt] = seq_wt
            out_dict[header_mut] = seq_mut

        else:
            logging.info(
                "[ESEfinder] Variant {} discarded. Reason: There is only support for SNVs".format(
                    id
                )
            )
        return out_dict

    def HEXplorer(
        self,
        record: cyvcf2.Variant,
        vep_annotation: str,
        variant_seq_map: dict,
        hexplorer_ref_seq: str,
        hexplorer_mut_seq: str,
    ):
        """
        121bp sequences will be generated with each
        variant being located in the middle position (61).
        """
        id = (
            record.CHROM
            + "_"
            + str(record.POS)
            + "_"
            + record.REF
            + "_"
            + record.ALT[0]
        )

        strand = vep_annotation.split("|")[self.vep_indexes["strand"]]
        if record.var_type in ["snp"]:
            if strand == "1":
                seq_wt = str(
                    self.fasta[record.CHROM][record.POS - 61 : record.POS + 60]
                )
                seq_mut = seq_wt[:60] + record.ALT[0] + seq_wt[61:]

            elif strand == "-1":
                seq_wt = str(
                    -self.fasta[record.CHROM][record.POS - 61 : record.POS + 60]
                )
                seq_mut = (
                    seq_wt[:60] + str(Seq(record.ALT[0]).complement()) + seq_wt[61:]
                )

            if len(variant_seq_map) >= 1:
                variant_seq_map[id] = 121 * len(variant_seq_map) + 121 - 61
            else:
                variant_seq_map[id] = 60

            hexplorer_ref_seq += seq_wt
            hexplorer_mut_seq += seq_mut

            return variant_seq_map, hexplorer_ref_seq, hexplorer_mut_seq

        else:
            logging.info(
                "[HEXplorer] Variant {} discarded. Reason: There is only support for SNVs".format(
                    id
                )
            )
            return variant_seq_map, hexplorer_ref_seq, hexplorer_mut_seq

    def SpliceRover(self, record: cyvcf2.Variant, vep_annotation: str, out_dict: dict):
        """
        400bp sequences will be generated with the putative
        splice sites around the middle positions (201,202)
        """
        assert self.ss in [
            "donor",
            "acceptor",
        ], "--splicerover requires --splice_site to be set (donor or acceptor)"
        id = (
            record.CHROM
            + "_"
            + str(record.POS)
            + "_"
            + record.REF
            + "_"
            + record.ALT[0]
        )
        header_wt = id + "_WT"
        header_mut = id + "_Mutated"
        strand = vep_annotation.split("|")[self.vep_indexes["strand"]]

        if record.var_type in ["snp"]:
            if strand == "1":
                seq_wt = str(self.fasta[record.CHROM][record.POS - 3 : record.POS + 1])
                seq_mut = seq_wt[:2] + record.ALT[0] + seq_wt[-1]

                if _no_ss(id, seq_wt, seq_mut, self.ss, "SpliceRover"):
                    return out_dict

                seq_wt = str(
                    self.fasta[record.CHROM][record.POS - 250 : record.POS + 250]
                )
                seq_mut = seq_wt[:249] + record.ALT[0] + seq_wt[250:]

            elif strand == "-1":
                # Just to check motif presence
                seq_wt = str(-self.fasta[record.CHROM][record.POS - 2 : record.POS + 2])
                seq_mut = (
                    seq_wt[:2] + Seq(record.ALT[0]).reverse_complement() + seq_wt[-1]
                )

                if _no_ss(id, seq_wt, seq_mut, self.ss, "SpliceRover"):
                    return out_dict

                seq_wt = str(
                    -self.fasta[record.CHROM][record.POS - 250 : record.POS + 250]
                )
                seq_mut = (
                    seq_wt[:249] + str(Seq(record.ALT[0]).complement()) + seq_wt[250:]
                )

            out_dict[header_wt] = seq_wt
            out_dict[header_mut] = seq_mut

        else:
            logging.info(
                "[SpliceRover {}] Variant {} discarded. Reason: There is only support for SNVs".format(
                    self.ss, id
                )
            )
        return out_dict

    def Spliceator(self, record: cyvcf2.Variant, vep_annotation: str, out_dict: dict):
        """
        200bp sequences will be generated with the putative
        splice sites around the middle positions (101,102).
        """

        assert self.ss in [
            "donor",
            "acceptor",
        ], "--spliceator requires --splice_site to be set (donor or acceptor)"
        id = (
            record.CHROM
            + "_"
            + str(record.POS)
            + "_"
            + record.REF
            + "_"
            + record.ALT[0]
        )
        header_wt = id + "_WT"
        header_mut = id + "_Mutated"
        strand = vep_annotation.split("|")[self.vep_indexes["strand"]]

        if record.var_type in ["snp"]:
            if strand == "1":
                seq_wt = str(self.fasta[record.CHROM][record.POS - 3 : record.POS + 1])
                seq_mut = seq_wt[:2] + record.ALT[0] + seq_wt[-1]

                if _no_ss(id, seq_wt, seq_mut, self.ss, "Spliceator"):
                    return out_dict

                seq_wt = str(
                    self.fasta[record.CHROM][record.POS - 100 : record.POS + 100]
                )
                seq_mut = seq_wt[:99] + record.ALT[0] + seq_wt[100:]

            elif strand == "-1":
                # Just to check motif presence
                seq_wt = str(-self.fasta[record.CHROM][record.POS - 2 : record.POS + 2])
                seq_mut = (
                    seq_wt[:2] + Seq(record.ALT[0]).reverse_complement() + seq_wt[-1]
                )

                if _no_ss(id, seq_wt, seq_mut, self.ss, "Spliceator"):
                    return out_dict

                seq_wt = str(
                    -self.fasta[record.CHROM][record.POS - 100 : record.POS + 100]
                )
                seq_mut = (
                    seq_wt[:99] + str(Seq(record.ALT[0]).complement()) + seq_wt[100:]
                )

            out_dict[header_wt] = seq_wt
            out_dict[header_mut] = seq_mut

        else:
            logging.info(
                "[Spliceator {}] Variant {} discarded. Reason: There is only support for SNVs".format(
                    self.ss, id
                )
            )
        return out_dict

    def DSSP(self, record: cyvcf2.Variant, vep_annotation: str, out_dict: dict):
        """
        140bp sequences will be generated with the
        putative splice sites at position 69/70 (for acceptors)
        and 71/72 (for donors)
        """
        assert self.ss in [
            "donor",
            "acceptor",
        ], "--dssp requires --splice_site to be set (donor or acceptor)"
        id = (
            record.CHROM
            + "_"
            + str(record.POS)
            + "_"
            + record.REF
            + "_"
            + record.ALT[0]
        )
        header_wt = id + "_WT"
        header_mut = id + "_Mutated"
        strand = vep_annotation.split("|")[self.vep_indexes["strand"]]
        hgvs = vep_annotation.split("|")[self.vep_indexes["hgvsc"]]

        if record.var_type in ["snp"]:
            if self.ss == "donor":
                ss = "GT"
                pos_in_seq_bef = 71
                pos_in_seq_after = 69
            elif self.ss == "acceptor":
                ss = "AG"
                pos_in_seq_bef = 69
                pos_in_seq_after = 71

            if strand == "1":
                seq_wt = str(self.fasta[record.CHROM][record.POS - 3 : record.POS + 1])
                seq_mut = seq_wt[:2] + record.ALT[0] + seq_wt[-1]

                if _no_ss(id, seq_wt, seq_mut, self.ss, "DSSP"):
                    return out_dict

                index = seq_mut.index(ss) if ss in seq_mut else seq_wt.index(ss)
                d = {0: 2, 1: 1, 2: 0}

                ss_pos = record.POS - d[index]
                mut_pos = pos_in_seq_bef + d[index]

                seq = str(
                    self.fasta[record.CHROM][
                        ss_pos - pos_in_seq_bef : ss_pos + pos_in_seq_after
                    ]
                )
                mut_seq = seq[: mut_pos - 1] + record.ALT[0] + seq[mut_pos:]

            elif strand == "-1":
                seq_wt = str(-self.fasta[record.CHROM][record.POS - 2 : record.POS + 2])
                seq_mut = (
                    seq_wt[:2] + Seq(record.ALT[0]).reverse_complement() + seq_wt[-1]
                )

                if _no_ss(id, seq_wt, seq_mut, self.ss, "DSSP"):
                    return out_dict

                index = seq_mut.index(ss) if ss in seq_mut else seq_wt.index(ss)
                d = {0: 1, 1: 0, 2: -1}
                ss_pos = record.POS + d[index]
                mut_pos = pos_in_seq_bef + d[index]
                seq = str(
                    -self.fasta[record.CHROM][
                        ss_pos - pos_in_seq_after : ss_pos + pos_in_seq_bef
                    ]
                )

                mut_seq = (
                    seq[:mut_pos]
                    + str(Seq(record.ALT[0]).complement())
                    + seq[mut_pos + 1 :]
                )

            if self.ss == "donor":
                assert any(
                    x[70:72] for x in [seq, mut_seq]
                ), "Positions 71 and 72 should contain GT"
            else:
                assert any(
                    x[68:70] for x in [seq, mut_seq]
                ), "Positions 69 and 70 should contain AG"

            out_dict[header_wt] = seq
            out_dict[header_mut] = mut_seq

        else:
            logging.info(
                "[Splice2Deep {}] Variant {} discarded. Reason: There is only support for SNVs".format(
                    self.ss, id
                )
            )
        return out_dict

    def Splice2Deep(self, record: cyvcf2.Variant, vep_annotation: str, out_dict: dict):
        assert self.ss in [
            "donor",
            "acceptor",
        ], "--splice2deep requires --ss to be set (donor or acceptor)"

        id = (
            record.CHROM
            + "_"
            + str(record.POS)
            + "_"
            + record.REF
            + "_"
            + record.ALT[0]
        )
        header_wt = id + "_WT"
        header_mut = id + "_Mutated"
        strand = vep_annotation.split("|")[self.vep_indexes["strand"]]
        hgvs = vep_annotation.split("|")[self.vep_indexes["hgvsc"]]

        if record.var_type in ["snp"]:
            if strand == "1":
                seq_wt = str(self.fasta[record.CHROM][record.POS - 3 : record.POS + 1])
                seq_mut = seq_wt[:2] + record.ALT[0] + seq_wt[-1]

                if _no_ss(id, seq_wt, seq_mut, self.ss, "Splice2Deep"):
                    return out_dict

                if self.ss == "donor":
                    ss = "GT"
                elif self.ss == "acceptor":
                    ss = "AG"

                index = seq_mut.index(ss) if ss in seq_mut else seq_wt.index(ss)
                d = {0: 2, 1: 1, 2: 0}

                ss_pos = record.POS - d[index]
                mut_pos = 300 + d[index]

                seq = str(self.fasta[record.CHROM][ss_pos - 1 - 300 : ss_pos + 1 + 300])
                mut_seq = seq[:mut_pos] + record.ALT[0] + seq[mut_pos + 1 :]

            elif strand == "-1":
                seq_wt = str(-self.fasta[record.CHROM][record.POS - 2 : record.POS + 2])
                seq_mut = (
                    seq_wt[:2] + Seq(record.ALT[0]).reverse_complement() + seq_wt[-1]
                )

                if _no_ss(id, seq_wt, seq_mut, self.ss, "Splice2Deep"):
                    return out_dict

                if self.ss == "donor":
                    ss = "GT"
                elif self.ss == "acceptor":
                    ss = "AG"

                # # index = 0: mutation is at the pos + 3 in the intron
                # # index = 1: mutation triggers the creation/lost of the T
                # # index = 2: mutation triggers the creation/lost of the G
                index = seq_mut.index(ss) if ss in seq_mut else seq_wt.index(ss)
                d = {0: 1, 1: 0, 2: -1}
                ss_pos = record.POS + d[index]
                mut_pos = 301 + d[index]
                seq = str(
                    -self.fasta[record.CHROM][ss_pos - 1 - 300 : ss_pos + 1 + 300]
                )
                mut_seq = (
                    seq[:mut_pos]
                    + str(Seq(record.ALT[0]).complement())
                    + seq[mut_pos + 1 :]
                )

            ss = "GT" if self.ss == "donor" else "AG"
            assert any(
                x[300:302] == ss for x in [seq, mut_seq]
            ), "Positions 301 and 302 should be {}".format(ss)
            out_dict[header_wt] = seq
            out_dict[header_mut] = mut_seq

        else:
            logging.info(
                "[Splice2Deep {}] Variant {} discarded. Reason: There is only support for SNVs".format(
                    self.ss, id
                )
            )
        return out_dict

    def MaxEntScan(self, record: cyvcf2.Variant, vep_annotation: str, out_dict: dict):
        """
        17bp (for donors) or 45bp (for acceptors) sequences will be generated
        to test the sliding window approach
        """
        assert self.ss in [
            "donor",
            "acceptor",
        ], "--maxentscan requires --splice_site to be set (donor or acceptor)"

        id = (
            record.CHROM
            + "_"
            + str(record.POS)
            + "_"
            + record.REF
            + "_"
            + record.ALT[0]
        )
        header_wt = id + "_WT"
        header_mut = id + "_Mutated"
        strand = vep_annotation.split("|")[self.vep_indexes["strand"]]
        context = 8 if self.ss == "donor" else 22
        if record.var_type in ["snp"]:
            if strand == "1":
                seq_wt = str(
                    self.fasta[record.CHROM][
                        record.POS - context - 1 : record.POS + context
                    ]
                )
                seq_mut = seq_wt[:context] + record.ALT[0] + seq_wt[context + 1 :]

            elif strand == "-1":
                seq_wt = str(
                    -self.fasta[record.CHROM][
                        record.POS - context - 1 : record.POS + context
                    ]
                )
                seq_mut = (
                    seq_wt[:context]
                    + str(Seq(record.ALT[0]).complement())
                    + seq_wt[context + 1 :]
                )
            out_dict[header_wt] = seq_wt
            out_dict[header_mut] = seq_mut

        else:
            logging.info(
                "[MaxEntScan {}] Variant {} discarded. Reason: There is only support for SNVs".format(
                    self.ss, id
                )
            )
        return out_dict

    def BranchpointTools(
        self,
        record: cyvcf2.Variant,
        vep_annotation: str,
        out_dict: dict,
        tool: str,
        seq_size: str = 500,
    ):
        hgvs = vep_annotation.split("|")[self.vep_indexes["hgvsc"]]

        v = self.hp.parse_hgvs_variant(hgvs.split(" ")[0])

        id = (
            record.CHROM
            + "_"
            + str(record.POS)
            + "_"
            + record.REF
            + "_"
            + record.ALT[0]
        )
        offset = v.posedit.pos.start.offset

        if offset == 0:
            logging.info(
                "[{}] Variant {} discarded. Reason: Not intronic ({}).".format(
                    tool, id, hgvs
                )
            )

        elif offset > 0:
            logging.info(
                "[{}] Variant {} discarded. Reason: Positive HGVSc offset ({}).".format(
                    tool, id, hgvs
                )
            )

        elif offset < -seq_size:
            logging.info(
                "[{}] Variant {} discarded. Reason: Too far away from splicing acceptor (Max={}; Observed={}).".format(
                    tool, id, seq_size, offset
                )
            )

        else:
            header_wt = id + "_WT"
            header_mut = id + "_Mutated"
            strand = vep_annotation.split("|")[self.vep_indexes["strand"]]

            if record.var_type in ["snp", "mnp"]:
                if strand == "1":
                    adjust_offset = 0 if record.var_type == "mnp" else -1
                    seq_down_wt = str(
                        self.fasta[record.CHROM][
                            record.POS - 1 : record.POS + abs(offset) + adjust_offset
                        ]
                    )
                    seq_down_mut = record.ALT[0] + seq_down_wt[len(record.REF) :]

                    bp_ups_wt = seq_size - len(seq_down_wt)
                    seq_upst = str(
                        self.fasta[record.CHROM][
                            record.POS - 1 - bp_ups_wt : record.POS - 1
                        ]
                    )

                elif strand == "-1":
                    seq_down_wt = str(
                        -self.fasta[record.CHROM][
                            record.POS
                            - abs(offset)
                            + len(record.REF)
                            - 1 : record.POS
                            + len(record.REF)
                            - 1
                        ]
                    )
                    seq_down_mut = (
                        str(Seq(record.ALT[0]).reverse_complement())
                        + seq_down_wt[len(record.REF) :]
                    )

                    bp_ups_wt = seq_size - len(seq_down_wt)
                    seq_upst = str(
                        -self.fasta[record.CHROM][
                            record.POS
                            + len(record.REF)
                            - 1 : record.POS
                            + bp_ups_wt
                            + len(record.REF)
                            - 1
                        ]
                    )

                out_dict[header_wt] = seq_upst + seq_down_wt
                out_dict[header_mut] = seq_upst + seq_down_mut

            elif record.var_subtype == "del":
                del_size = len(record.REF) - len(record.ALT[0])
                if strand == "1":
                    seq_down_wt = str(
                        self.fasta[record.CHROM][record.POS : record.POS + abs(offset)]
                    )
                    seq_down_mut = seq_down_wt[del_size:]

                    bp_ups_wt = seq_size - len(seq_down_wt)
                    bp_ups_mut = seq_size - len(seq_down_mut)

                    seq_upst_wt = str(
                        self.fasta[record.CHROM][
                            record.POS - 1 - bp_ups_wt : record.POS - 1
                        ]
                    )
                    seq_upst_mut = str(
                        self.fasta[record.CHROM][
                            record.POS - 1 - bp_ups_mut : record.POS - 1
                        ]
                    )

                elif strand == "-1":
                    seq_down_wt = str(
                        -self.fasta[record.CHROM][
                            record.POS - abs(offset) + del_size : record.POS + del_size
                        ]
                    )
                    seq_down_mut = seq_down_wt[del_size:]

                    bp_ups_wt = seq_size - len(seq_down_wt)
                    bp_ups_mut = seq_size - len(seq_down_mut)

                    seq_upst_wt = str(
                        -self.fasta[record.CHROM][
                            record.POS + del_size : record.POS + del_size + bp_ups_wt
                        ]
                    )
                    seq_upst_mut = str(
                        -self.fasta[record.CHROM][
                            record.POS + del_size : record.POS + del_size + bp_ups_mut
                        ]
                    )

                out_dict[header_wt] = seq_upst_wt + seq_down_wt
                out_dict[header_mut] = seq_upst_mut + seq_down_mut

            else:
                raise ValueError("Not seen case")

        return out_dict


def main():
    """
    Main function
    """
    parser = argparse.ArgumentParser(
        prog="vcf2seq",
        description="Tool to generate proper input for different set of splicing-related models",
    )

    parser.add_argument(
        dest="vcf",
        help="Input VCF file annotated with VEP, such that the HGVSc field exists in the output.",
    )
    parser.add_argument(
        dest="fasta",
        help="Reference genome in fasta format. Should be the same genome assembly as the variants represented in the vcf",
    )
    parser.add_argument(dest="outbasename", help="Outbasename to write the output")
    parser.add_argument(
        "--ss",
        type=str,
        choices=("donor", "acceptor"),
        help="For splice site prediction tools, refer to which splice site type VCF refers to",
    )
    parser.add_argument(
        "--vep_tag",
        default="CSQ",
        choices=["ANN", "CSQ"],
        help='Field in VCF where VEP annotations are stored. Default: "CSQ".',
    )
    parser.add_argument(
        "--surrounding",
        type=int,
        metavar="",
        help="Instead of extracting input for a specific tool, this argument allows extracting N number of base-pairs on each side of the variant.",
    )
    parser.add_argument(
        "--svm_bp_finder",
        action="store_true",
        help="Generate input for SVM-BP finder tool (500bp upstream of acceptor site). Only sequences that refer to the end of introns will be written",
    )
    parser.add_argument(
        "--bpp",
        action="store_true",
        help="Generate input for BPP (200bp upstream of acceptor site). Only sequences that refer to the end of introns will be written",
    )
    parser.add_argument(
        "--maxentscan",
        action="store_true",
        help="Generate input for maxentscanpy utility.",
    )
    parser.add_argument(
        "--splice2deep",
        action="store_true",
        help="Generate input for Splice2Deep donor model (602bp sequences). Putative splice site positions need to be in positions 300 and 301.",
    )
    parser.add_argument(
        "--spliceator",
        action="store_true",
        help="Generate input for Spliceator acceptor model (200bp sequences). Putative splice site positions will be around positions 100 and 101.",
    )
    parser.add_argument(
        "--dssp",
        action="store_true",
        help="Generate input for DSSP model (140bp sequences). Putative splice site positions need be in positions 69 and 70 (for acceptors) and 71 and 72 (for donors).",
    )
    parser.add_argument(
        "--splicerover",
        action="store_true",
        help="Generate input for SpliceRover model (400bp sequences). Putative splice site positions will be around positions 200 and 201. Be default, we generate 500bp sequence so that splice sites can be predicted between position 150 and 250",
    )
    parser.add_argument(
        "--hexplorer",
        action="store_true",
        help="Generate input for HEXplorer web tool. Sequences of 121bp for each mutation will be generated (mutation at pos 61).",
    )
    parser.add_argument(
        "--esefinder",
        action="store_true",
        help="Generate input for ESEfinder web tool. Sequences of 21bp for each mutation will be generated (mutation at pos 11).",
    )
    parser.add_argument(
        "--esrseq",
        action="store_true",
        help="Generate input to calculate delta ESRseq scores. Sequences of 11bp for each mutation will be generated (mutation at pos 6)",
    )
    args = parser.parse_args()

    # Initial checks
    vep_indexes = validate_vep(args.vcf)

    # Iterate over VCF
    dp = DataProcessing(
        args.vcf,
        vep_indexes,
        args.vep_tag,
        fasta=args.fasta,
        outbasename=args.outbasename,
        ss=args.ss,
        svm_bp_finder=args.svm_bp_finder,
        bpp=args.bpp,
        maxentscan=args.maxentscan,
        splice2deep=args.splice2deep,
        spliceator=args.spliceator,
        dssp=args.dssp,
        splicerover=args.splicerover,
        hexplorer=args.hexplorer,
        esefinder=args.esefinder,
        esrseq=args.esrseq,
        surrounding=args.surrounding,
    )


if __name__ == "__main__":
    main()
