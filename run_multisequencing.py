import argparse

from multisequencing import MultiAligner

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--seqA_path", type=str, required=True,
                        help="Path to .fasta file.")
    parser.add_argument("--seqB_path", type=str, required=True,
                        help="Path to .fasta file.")
    parser.add_argument("--aligner_mode", type=str, default='global',
                        help="local or global.")
    parser.add_argument("--distance_matrix", type=str, default='distance_matrix.csv',
                        help="Distance matrix path.")
    args = parser.parse_args()

    multialigner = MultiAligner([args.seqA_path, args.seqB_path])

    multialigner.load_files()

    multialigner.get_alignments()

    multialigner.get_profile_matrix()

    multialigner.get_profile_alignment()

    print("debug")

