import argparse

from multisequencing import MultiAligner

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--seqA_path", type=str, required=True,
                        help="Path to first .fasta file.")
    parser.add_argument("--seqB_path", type=str, required=True,
                        help="Path to second .fasta file.")
    parser.add_argument("--distance_matrix", type=str, default='./clustal_matrix',
                        help="Distance matrix path.")
    parser.add_argument("--gap_penalty", type=bool, default=True,
                        help="True if gap penalty enabled.")
    parser.add_argument("--custom_matrix", type=bool, default=True,
                        help="True if want to use custom matrix")
    args = parser.parse_args()

    multialigner = MultiAligner([args.seqA_path, args.seqB_path], args.gap_penalty, args.custom_matrix,
                                args.distance_matrix)

    multialigner.load_files()

    multialigner.get_alignments()

    multialigner.get_profile_matrix()

    multialigner.get_slow_alignment()

    multialigner.get_progressive_profile_alignment()
