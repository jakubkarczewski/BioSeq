from os.path import isfile
from pickle import load

import pandas as pd
from Bio import Align, pairwise2


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
        self.aligner_mode = aligner_mode
        # self.aligner_open_gap_score = -1
        # self.aligner_extend_gap_score = 0

        # set flags for sequence format
        self.rna = rna
        self.amino = amino

        # for aminoacid_2_dna translation
        self.performed_transformations = {}

    def read_seqs(self, fasta=False):
        """Read sequences from files of specified format."""
        if not fasta:
            for path, name in zip(self.seq_paths, ('seqA', 'seqB')):
                with open(path) as f:
                    sequence = f.read().replace('\n', '') if not self.rna else self._rna_2_dna(f.read().
                                                                                               replace('\n', ''))
                    self.seqs[name] = sequence if not self.amino else self._dna_2_aminoacid(sequence)
        else:
            self._read_fasta()

    def get_score(self):
        """Returns best alignment score"""
        return self.aligner.score(**self.seqs)

    def _get_alignments(self, custom_matrix=False):
        """Returns alignments"""
        if custom_matrix:
            self._load_distance_matrix()

        kwargs = dict(self.seqs)
        kwargs['one_alignment_only'] = True
        kwargs['penalize_extend_when_opening'] = True
        kwargs['matrix'] = self.distance_matrix_dict

        if self.aligner_mode == 'local':
            return pairwise2.align.localds(kwargs['seqA'], kwargs['seqB'], self.distance_matrix_dict, -10, 0)
        elif self.aligner_mode == 'global':
            pass
        else:
            raise Exception("wrong mode.")

    def print_alignments(self,custom_matrix=False):
        """Print alignment and it's score"""

        for elem in self._get_alignments(custom_matrix=custom_matrix):
            print(elem)
            print(elem.___dict__)

    # def _get_best_alignment(self, custom_matrix=False):
    #     """Returns best alignment"""
    #     if custom_matrix:
    #         self._set_distance_matrix()
    #     kwargs = dict(self.seqs)
    #     kwargs['score_only'] = True
    #     return self.aligner.align(**kwargs)

    @classmethod
    def _rna_2_dna(cls, sequence):
        """Transforms RNA sequence do DNA sequence."""
        return sequence.replace('U', 'T')

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

    def _aminoacid_2_dna(self, sequence):
        raise NotImplementedError

    def _matrix_2_sequence(self, matrix):
        raise NotImplementedError

    def _get_distance(self, score):
        raise NotImplementedError

    def _read_fasta(self):
        raise NotImplementedError

    def _load_distance_matrix(self):
        """Load distance matrix dict from .csv file."""
        raw_dict = pd.read_csv(self.distance_matrix_path, index_col=0).to_dict()
        self.distance_matrix_dict = {}
        for key in raw_dict.keys():
            for subkey in raw_dict[key].keys():
                if str(raw_dict[key][subkey]) != 'nan' and subkey != '-' and key != '-':
                    self.distance_matrix_dict[(key, subkey)] = float(raw_dict[key][subkey])
