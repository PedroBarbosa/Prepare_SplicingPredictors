from operator import is_
from pickletools import read_uint1
import sys
from vcf_utils import *
import cyvcf2
import hgvs.parser
from Bio.Seq import Seq
from pyfaidx import Fasta
import argparse
from curses import KEY_SAVE
from typing import Union, TextIO, Optional
import pathlib
import logging
from numpy import rec
logging.basicConfig(stream=sys.stdout, level=logging.INFO, format='%(asctime)s %(message)s')


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
    def __init__(self, vcf: str,
                 vep_indexes: dict,
                 vep_tag: str,
                 fasta: Fasta,
                 outbasename: str,
                 svm_bp_finder: bool = False,
                 splice2deep_donor: bool = False,
                 splice2deep_acceptor: bool = False):

        self.hp = hgvs.parser.Parser()
        self.vep_tag = vep_tag
        self.vep_indexes = vep_indexes
        self.fasta = Fasta(fasta)
        self.outbasename = outbasename
        self.svm_bp_finder = svm_bp_finder
        self.splice2deep_donor = splice2deep_donor
        self.splice2deep_acceptor = splice2deep_acceptor
        self._iterate_vcf(vcf)

    def _iterate_vcf(self, vcf):

        fasta_svm_bp, fasta_splice2deep_donor, fasta_splice2deep_acceptor = {}, {}, {}

        for record in VCF(vcf):

            if len(record.ALT) > 1:
                raise ValueError(
                    'Multi allelic records are not allowed. Please split them before.')
            vep_annotation = record.INFO.get(self.vep_tag)

            if vep_annotation:
                # First consequence is selected
                vep_annotation = vep_annotation.split(",")[0]
            else:
                continue

            if self.svm_bp_finder:
                fasta_svm_bp = self.SVM_BP_finder(record,
                                                  vep_annotation,
                                                  out_dict=fasta_svm_bp)

            if self.splice2deep_donor:
                fasta_splice2deep_donor = self.Splice2Deep(record,
                                                           vep_annotation,
                                                           out_dict=fasta_splice2deep_donor,
                                                           is_donor=True)
            if self.splice2deep_acceptor:
                fasta_splice2deep_acceptor = self.Splice2Deep(record,
                                                              vep_annotation,
                                                              out_dict=fasta_splice2deep_acceptor,
                                                              is_donor=False)

        if self.svm_bp_finder:
            _dict_to_fasta(
                fasta_svm_bp, '{}_SVM_BP_finder.fa'.format(self.outbasename))

        if self.splice2deep_donor:
            _dict_to_fasta(
                fasta_splice2deep_donor, '{}_Splice2Deep_donor.fa'.format(self.outbasename))
        
        if self.splice2deep_acceptor:
            _dict_to_fasta(
                fasta_splice2deep_acceptor, '{}_Splice2Deep_acceptor.fa'.format(self.outbasename))
            
    def Splice2Deep(self,
                    record: cyvcf2.Variant,
                    vep_annotation: str,
                    out_dict: dict,
                    is_donor: bool):

        def _no_donor(seq_wt, seq_mut):
            if all("GT" not in x for x in [seq_wt, seq_mut]):
                logging.info('[Splice2Deep donor] Variant {} discarded. Reason: no GT dinucleotide formed at the 4-mer around the mutation (wt={}, mut={}).'.format(id, seq_wt, seq_mut))    
                return True
            else:
                return False

        def _no_acceptor(seq_wt, seq_mut):
            if all("AG" not in x for x in [seq_wt, seq_mut]):
                logging.info('[Splice2Deep acceptor] Variant {} discarded. Reason: no AG dinucleotide formed at the 4-mer around the mutation (wt={}, mut={}).'.format(id, seq_wt, seq_mut))    
                return True
            else:
                return False
            
        id = record.CHROM + "_" + \
            str(record.POS) + "_" + record.REF + "_" + record.ALT[0]
        header_wt = id + "_WT"
        header_mut = id + "_Mutated"
        strand = vep_annotation.split("|")[self.vep_indexes['strand']]
        hgvs = vep_annotation.split("|")[self.vep_indexes['hgvsc']]
        
        if record.var_type in ['snp']:
            if strand == "1":
                seq_wt = str(self.fasta[record.CHROM][record.POS - 3:record.POS + 1])
                seq_mut = seq_wt[:2] + record.ALT[0] + seq_wt[-1]
                
                if is_donor:
                    if _no_donor(seq_wt, seq_mut):
                        return out_dict
                    ss="GT"
                else:
                    if _no_acceptor(seq_wt, seq_mut):
                        return out_dict
                    ss="AG"
                    
                index = seq_mut.index(ss) if ss in seq_mut else seq_wt.index(ss)
                d = {0: 2, 1: 1, 2: 0}

                ss_pos = record.POS - d[index]
                mut_pos = 300 + d[index]

                seq = str(self.fasta[record.CHROM][ss_pos -1 -300:ss_pos + 1 + 300])
                mut_seq = seq[:mut_pos] + record.ALT[0] + seq[mut_pos + 1:]
       
            elif strand == "-1":
                seq_wt = str(-self.fasta[record.CHROM][record.POS -2:record.POS + 2])
                seq_mut = seq_wt[:2] + Seq(record.ALT[0]).reverse_complement() + seq_wt[-1]

                if is_donor:
                    if _no_donor(seq_wt, seq_mut):
                        return out_dict
                    ss="GT"
                else:
                    if _no_acceptor(seq_wt, seq_mut):
                        return out_dict
                    ss="AG"
                    
                # # index = 0: mutation is at the pos + 3 in the intron
                # # index = 1: mutation triggers the creation/lost of the T 
                # # index = 2: mutation triggers the creation/lost of the G
                index = seq_mut.index(ss) if ss in seq_mut else seq_wt.index(ss)
                d = {0: 1, 1: 0, 2: -1}
                ss_pos = record.POS + d[index]
                mut_pos = 301 + d[index]
                seq = str(-self.fasta[record.CHROM][ss_pos -1 -300:ss_pos + 1 + 300])
                mut_seq = seq[:mut_pos] + str(Seq(record.ALT[0]).complement()) + seq[mut_pos + 1:]

            ss = "GT" if is_donor else "AG"
            assert any(x[300:302] == ss for x in [seq, mut_seq]), "Positions 301 and 302 should be {}".format(ss)
            out_dict[header_wt] = seq
            out_dict[header_mut] = mut_seq
           
        else:
            logging.info("[Splice2Deep {}] Variant {} discarded. Reason: There is only support for SNVs".format('donor' if is_donor else 'acceptor', id))
        return out_dict
     
    def SVM_BP_finder(self,
                      record: cyvcf2.Variant,
                      vep_annotation: str,
                      out_dict: dict):

        hgvs = vep_annotation.split("|")[self.vep_indexes['hgvsc']]

        v = self.hp.parse_hgvs_variant(hgvs.split(" ")[0])

        id = record.CHROM + "_" + \
            str(record.POS) + "_" + record.REF + "_" + record.ALT[0]
        offset = v.posedit.pos.start.offset

        if offset == 0:
            logging.info(
                '[SVM-BP_finder] Variant {} discarded. Reason: Not intronic ({}).'.format(id, hgvs))

        elif offset > 0:
            logging.info(
                '[SVM-BP_finder] Variant {} discarded. Reason: Positive HGVSc offset ({}).'.format(id, hgvs))

        elif offset < -500:
            logging.info(
                '[SVM-BP_finder] Variant {} discarded. Reason: Too far away from splicing acceptor (Max=500; Observed={}).'.format(id, offset))

        else:
            header_wt = id + "_WT"
            header_mut = id + "_Mutated"
            strand = vep_annotation.split("|")[self.vep_indexes['strand']]

            if record.var_type in ['snp', 'mnp']:

                if strand == "1":
                    adjust_offset = 0 if record.var_type == 'mnp' else - 1
                    seq_down_wt = str(
                        self.fasta[record.CHROM][record.POS - 1:record.POS + abs(offset) + adjust_offset])
                    seq_down_mut = record.ALT[0] + \
                        seq_down_wt[len(record.REF):]

                    bp_ups_wt = 500 - len(seq_down_wt)
                    seq_upst = str(
                        self.fasta[record.CHROM][record.POS - 1 - bp_ups_wt: record.POS - 1])

                elif strand == "-1":

                    seq_down_wt = str(-self.fasta[record.CHROM][record.POS - abs(
                        offset) + len(record.REF) - 1: record.POS + len(record.REF) - 1])
                    seq_down_mut = str(
                        Seq(record.ALT[0]).reverse_complement()) + seq_down_wt[len(record.REF):]

                    bp_ups_wt = 500 - len(seq_down_wt)
                    seq_upst = str(-self.fasta[record.CHROM][record.POS + len(
                        record.REF) - 1: record.POS + bp_ups_wt + len(record.REF) - 1])

                out_dict[header_wt] = seq_upst + seq_down_wt
                out_dict[header_mut] = seq_upst + seq_down_mut

            elif record.var_subtype == 'del':

                del_size = len(record.REF) - len(record.ALT[0])
                if strand == "1":

                    seq_down_wt = str(
                        self.fasta[record.CHROM][record.POS: record.POS + abs(offset)])
                    seq_down_mut = seq_down_wt[del_size:]

                    bp_ups_wt = 500 - len(seq_down_wt)
                    bp_ups_mut = 500 - len(seq_down_mut)

                    seq_upst_wt = str(
                        self.fasta[record.CHROM][record.POS - 1 - bp_ups_wt: record.POS - 1])
                    seq_upst_mut = str(
                        self.fasta[record.CHROM][record.POS - 1 - bp_ups_mut: record.POS - 1])

                elif strand == "-1":

                    seq_down_wt = str(-self.fasta[record.CHROM][record.POS - abs(
                        offset) + del_size: record.POS + del_size])
                    seq_down_mut = seq_down_wt[del_size:]

                    bp_ups_wt = 500 - len(seq_down_wt)
                    bp_ups_mut = 500 - len(seq_down_mut)

                    seq_upst_wt = str(-self.fasta[record.CHROM][record.POS +
                                      del_size: record.POS + del_size + bp_ups_wt])
                    seq_upst_mut = str(-self.fasta[record.CHROM][record.POS +
                                       del_size: record.POS + del_size + bp_ups_mut])

                out_dict[header_wt] = seq_upst_wt + seq_down_wt
                out_dict[header_mut] = seq_upst_mut + seq_down_mut

            else:
                raise ValueError('Not seen case')

        return out_dict


def main():
    """
    Main function
    """
    parser = argparse.ArgumentParser(prog="vcf2seq",
                                     description='Tool to generate proper input for different set of splicing-related models')

    parser.add_argument(
        dest='vcf', help='Input VCF file annotated with VEP, such that the HGVSc field exists in the output.')
    parser.add_argument(
        dest='fasta', help='Reference genome in fasta format. Should be the same genome assembly as the variants represented in the vcf')
    parser.add_argument(dest='outbasename',
                        help='Outbasename to write the output')

    parser.add_argument('--vep_tag', default='CSQ', choices=[
                        'ANN', 'CSQ'], help='Field in VCF where VEP annotations are stored. Default: "CSQ".')
    parser.add_argument('--svm_bp_finder', action='store_true',
                        help='Generate input for SVM-BP finder tool. Only sequences that refer to the end of introns will be written')
    parser.add_argument('--splice2deep_donor', action='store_true',
                        help='Generate input for Splice2Deep donor model (602bp sequences). Putative splice site positions should be in positions 300 and 301.')
    parser.add_argument('--splice2deep_acceptor', action='store_true',
                        help='Generate input for Splice2Deep acceptor model (602bp sequences). Putative splice site positions should be in positions 300 and 301.')
    args = parser.parse_args()

    # Initial checks
    vep_indexes = validate_vep(args.vcf)

    # Iterate over VCF
    dp = DataProcessing(args.vcf,
                        vep_indexes,
                        args.vep_tag,
                        fasta=args.fasta,
                        outbasename=args.outbasename,
                        svm_bp_finder=args.svm_bp_finder,
                        splice2deep_donor=args.splice2deep_donor,
                        splice2deep_acceptor=args.splice2deep_acceptor)


if __name__ == '__main__':
    main()
