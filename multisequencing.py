"""
Zestawienia wielu sekwencji. Typ A.
Dla zadanego wielodopasowania ciągów
nukleotydowych i macierzy podobieństwa wyznaczyć macierz profilu i słowo konsensusowe.

Dla dwóch wielodopasowań:
- zestawić je w jedno przez złożenie (optymalne, globalne, bez
funkcji kary za przerwy, macierz podobieństwa liter jest dana) ich profili.

Wreszcie:
 - progressive multialigning podanych sekwencji w jedno wielodopasowanie poprzez
klasteryzację i zestawienie profili. Strategia wyboru sklejanych pod-dopasowań: dowolna
prosta np. UPGMA (tj. guide tree nie wymagane).
"""
import os
from os.path import join, isfile

import numpy as np
from Bio import AlignIO
from Bio.Align.Applications import ClustalwCommandline
from Bio.Alphabet import generic_dna
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from Bio.Align import MultipleSeqAlignment, AlignInfo

from sequencing import Aligner

class MultiAligner:
    def __init__(self, seq_path, clustal_path="./clustalw2"):
        assert isfile(seq_path)
        assert  isfile(clustal_path)
        self.clustal_path = clustal_path
        self.seqs_path = seq_path
        self.seqs = None
        self.alignments = None

    def open_phy_file(self):
        """Opens input file and loads it's content."""
        self.seqs = AlignIO.parse(self.seqs_path, "fasta")

    def get_alignments(self):
        cline = ClustalwCommandline(self.clustal_path, infile=self.seqs_path)
        stdout, stderr = cline()
        if not stderr:
            print(stdout)
        else:
            print(stderr)

    @staticmethod
    def normalize(one_dict):
        rv = {}
        one_sum = sum(list(one_dict.values()))
        for key in one_dict.keys():
            rv[key] = one_dict[key] / one_sum
        return rv

    def process_matrix(self, pssm):
        processed_list = []
        for i, elem in enumerate(pssm):
            processed_list.append((str(i), self.normalize(elem[-1])))
        return  processed_list

    def get_profile_matrix(self):
        filename = self.seqs_path.split('.')[0] + ".aln"
        alignment = AlignIO.read(filename, 'clustal')
        summary_align = AlignInfo.SummaryInfo(alignment)
        consensus = summary_align.dumb_consensus()
        my_pssm = summary_align.pos_specific_score_matrix(consensus)
        copy = my_pssm.pssm
        processed = self.process_matrix(copy)
        my_pssm.pssm = processed
        print(my_pssm)
        print(consensus)


    # @staticmethod
    # def _preprocess(seqs):
    #     preprocessed = []
    #     for seq in seqs:
    #         preprocessed.append(SeqRecord(Seq(seq, generic_dna)))
    #     return preprocessed
    #     # return [SeqRecord(Seq(seq, generic_dna)) for seq in seqs]
    #
    # def multiple_alignment(self):
    #     self.alignments = MultipleSeqAlignment(self._preprocess(self.seqs))
