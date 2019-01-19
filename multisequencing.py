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
    def __init__(self, seq_paths, clustal_path="./clustalw2"):
        for path in seq_paths:
            assert isfile(path)
        assert isfile(clustal_path)
        self.clustal_path = clustal_path
        self.seqs_paths = seq_paths
        self.seqs = {'seqA': None, 'seqB': None}
        self.alignments = {'seqA': None, 'seqB': None}

    def load_files(self):
        """Opens input file and loads it's content."""
        for path in self.seqs_paths:
            seq_name = path.split('.')[0]
            self.seqs[seq_name] = AlignIO.parse(path, "fasta")

    def get_alignments(self, verbose=False):
        for path in self.seqs_paths:
            cline = ClustalwCommandline(self.clustal_path, infile=path)
            stdout, stderr = cline()
            if verbose:
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
        return processed_list

    def get_profile_matrix(self):
        for path in self.seqs_paths:
            filename = path.split('.')[0] + ".aln"
            alignment = AlignIO.read(filename, 'clustal')
            summary_align = AlignInfo.SummaryInfo(alignment)
            consensus = summary_align.dumb_consensus()
            my_pssm = summary_align.pos_specific_score_matrix(consensus)
            copy = my_pssm.pssm
            processed = self.process_matrix(copy)
            my_pssm.pssm = processed
            print('For', path.split('.')[0], ':')
            print("Profile matrix:\n", my_pssm)
            print("Consensus word:\n", consensus)

    def get_profile_alignment(self):
        filenames = [base + '.aln' for base in self.seqs.keys()]
        assert len(filenames) == 2
        cline = ClustalwCommandline(self.clustal_path, profile1=filenames[0], profile2=filenames[1])
        stdout, stderr = cline()

        print(stderr, '\n', stdout)




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
