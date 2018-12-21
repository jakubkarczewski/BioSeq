from os.path import isfile
from pickle import load

import pandas as pd
from Bio import Align
from Levenshtein import distance


class Aligner:
    def __init__(self, seq_paths, aminoacid_dict_path, aligner_mode, distance_matrix_path, rna=False, amino=False):
        # check validity of data
        assert isinstance(seq_paths, list)
        assert all(isfile(path) for path in seq_paths)
        assert len(seq_paths) == 2
        assert aligner_mode in ('local', 'global')
        assert isfile(aminoacid_dict_path)
        assert isfile(distance_matrix_path)

        self.seq_paths = seq_paths
        self.seqs = {}
        self.distance_matrix_path = distance_matrix_path
        self.distance_matrix_dict = None

        # load dna_2_aminoacid dict
        with open(aminoacid_dict_path, 'rb') as f:
            self.aminoacid_dict = load(f)

        # create aligner
        self.aligner = Align.PairwiseAligner()
        self.aligner.mode = aligner_mode
        self.aligner.open_gap_score = -1
        self.aligner.extend_gap_score = 0

        # set flags for sequence format
        self.rna = rna
        self.amino = amino

        # for aminoacid_2_dna translation
        self.performed_transformations = {}

    def read_seqs(self):
        """Read sequences from files of specified format."""
        for path, name in zip(self.seq_paths, ('seqA', 'seqB')):
            with open(path) as f:
                sequence = f.read().replace('\n', '')
                self.seqs[name] = sequence

    def get_score(self):
        """Returns best alignment score"""
        return self.aligner.score(**self.seqs)

    def get_alignments(self, custom_matrix=False, analyse_amino=False):
        """Returns alignments"""
        if custom_matrix:
            self._load_distance_matrix()
            self.aligner.substitution_matrix = self.distance_matrix_dict

        if analyse_amino:
            return self.aligner.align(self._dna_2_aminoacid(self.seqs['seqA']),
                                      self._dna_2_aminoacid(self.seqs['seqB']))

        return self.aligner.align(**self.seqs)

    def rna_2_dna(self):
        """Transforms RNA sequence do DNA sequence."""
        temp_seqs = dict(self.seqs)
        for seq_name in self.seqs.keys():
            self.seqs[seq_name] = temp_seqs[seq_name].replace('U', 'T')

    def _dna_2_rna(self):
        """Transforms DNA sequence do RNA sequence."""
        temp_seqs = dict(self.seqs)
        for seq_name in self.seqs.keys():
            self.seqs[seq_name] = temp_seqs[seq_name].replace('T', 'U')

    def _dna_2_aminoacid(self, sequence):
        """Transforms DNA sequence to aminoacid sequence."""

        if len(sequence) % 3 != 0:
            sequence = sequence[:-(len(sequence) % 3)]

        translated = []
        for i in range(0, len(sequence), 3):
            trio = sequence[i:i + 3]
            aminoacid = self.aminoacid_dict[trio]
            translated.append(aminoacid)
            self.performed_transformations[aminoacid] = trio

        return ''.join(translated)

    def _aminoacid_2_dna(self, amino_sequence):
        """Transforms aminoacid sequence to DNA sequence."""
        return ''.join([self.performed_transformations[elem] for elem in amino_sequence])

    def _load_distance_matrix(self):
        """Load distance matrix dict from .csv file."""
        raw_dict = pd.read_csv(self.distance_matrix_path, index_col=0).to_dict()
        self.distance_matrix_dict = {}
        for key in raw_dict.keys():
            for subkey in raw_dict[key].keys():
                if str(raw_dict[key][subkey]) != 'nan' and subkey != '-' and key != '-':
                    self.distance_matrix_dict[(key, subkey)] = float(raw_dict[key][subkey])

    def get_levenstein_distance(self, analyse_amino=False):
        """Returns Levenstein edit distance"""
        return distance(*self.seqs.values()) if not analyse_amino else distance(self._dna_2_aminoacid(self.seqs['seqA']),
                                                                        self._dna_2_aminoacid(self.seqs['seqB']))
