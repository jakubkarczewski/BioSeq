from os.path import isfile
from pickle import load

from Bio import Align


class Aligner:
    def __init__(self, seq_paths, aminoacid_dict_path, aligner_mode, rna=False, amino=False):
        # check validity of data
        assert isinstance(seq_paths, list)
        assert all(isfile(path) for path in seq_paths)
        assert len(seq_paths) == 2
        assert aligner_mode in ('local', 'global')
        assert isfile(aminoacid_dict_path)

        self.seq_paths = seq_paths
        self.seqs = {}

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

    def _get_alignments(self):
        """Returns alignments"""
        return self.aligner.align(**self.seqs)

    def print_alignment(self):
        """Print alignment and it's score"""
        for elem in self._get_alignments():
            print(elem)
            print(elem.score)

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

