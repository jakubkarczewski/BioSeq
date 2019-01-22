from os.path import isfile

from Bio import AlignIO
from Bio.Align.Applications import ClustalwCommandline
from Bio.Align import AlignInfo


class MultiAligner:
    def __init__(self, seq_paths, gap_penalty, custom_matrix, matrix_path,  clustal_path="./clustalw2"):
        for path in seq_paths:
            assert isfile(path)
        assert isfile(clustal_path)
        self.clustal_path = clustal_path
        self.seqs_paths = seq_paths

        self.seqs = {'seqA': None, 'seqB': None}
        self.alignments = {'seqA': None, 'seqB': None}
        self.gap_penalty = gap_penalty
        self.custom_matrix = custom_matrix

        assert isfile(matrix_path)
        self.matrix_path = matrix_path

    def load_files(self):
        """Opens input files and loads their content."""
        for path in self.seqs_paths:
            seq_name = path.split('.')[0]
            self.seqs[seq_name] = AlignIO.parse(path, "fasta")

    def get_alignments(self):
        """Performs multialignment for sequences in each input file."""
        for path in self.seqs_paths:
            # set command parameters for clustalw2 execution
            args = [self.clustal_path]
            kwargs = {'infile': path}
            if self.custom_matrix:
                kwargs['transweight'] = 0
                kwargs['matrix'] = self.matrix_path
            # execute command
            cline = ClustalwCommandline(*args, **kwargs)
            stdout, stderr = cline()

            print(stdout, '\n', stderr)

    @staticmethod
    def normalize(one_dict):
        rv = {}
        one_sum = sum(list(one_dict.values()))
        for key in one_dict.keys():
            rv[key] = one_dict[key] / one_sum
        return rv

    def process_matrix(self, pssm):
        """Normalizes PSSM matrix."""
        processed_list = []
        for i, elem in enumerate(pssm):
            processed_list.append((str(i), self.normalize(elem[-1])))
        return processed_list

    def get_profile_matrix(self):
        """Prints PSSM matrix, consensus word."""
        for path in self.seqs_paths:
            filename = path.split('.')[0] + ".aln"
            alignment = AlignIO.read(filename, 'clustal')
            # get consensus word
            summary_align = AlignInfo.SummaryInfo(alignment)
            consensus = summary_align.dumb_consensus()
            # get pssm matrix and process it
            my_pssm = summary_align.pos_specific_score_matrix(consensus)
            copy = my_pssm.pssm
            processed = self.process_matrix(copy)
            my_pssm.pssm = processed
            # print the data
            print('For', path.split('.')[0], ':')
            print("Profile matrix:\n", my_pssm)
            print("Consensus word:\n", consensus)

    def get_slow_alignment(self):
        """Performs slow (non-progressive) multialignment."""
        filenames = [base + '.aln' for base in self.seqs.keys()]
        assert len(filenames) == 2
        # set command parameters for clustalw2 execution
        args = [self.clustal_path]
        kwargs = {'profile1': filenames[0], 'profile2': filenames[1]}
        kwargs['pwgapopen'] = 5 if self.gap_penalty else 0
        kwargs['pwgapext'] = 5 if self.gap_penalty else 0
        if self.custom_matrix:
            kwargs['pwmatrix'] = self.matrix_path
        # execute command
        cline = ClustalwCommandline(*args, **kwargs)
        stdout, stderr = cline()

        print(stdout, '\n', stderr)

    def get_progressive_profile_alignment(self):
        args = [self.clustal_path]
        kwargs = {'profile1': self.seqs_paths[0], 'profile2': self.seqs_paths[1]}
        kwargs['gapopen'] = 5 if self.gap_penalty else 0
        kwargs['gapext'] = 5 if self.gap_penalty else 0
        if self.custom_matrix:
            kwargs['transweight'] = 0
            kwargs['matrix'] = self.matrix_path

        cline = ClustalwCommandline(*args, **kwargs)
        stdout, stderr = cline()

        print(stdout, '\n', stderr)
