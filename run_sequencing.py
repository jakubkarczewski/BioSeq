import argparse

from sequencing import Aligner

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--seq_A_path", type=str, required=True,
                        help="Path to first seq file.")
    parser.add_argument("--seq_B_path", type=str, required=True,
                        help="Path to second seq file.")
    parser.add_argument("--aminoacid_dict_path", type=str, required=True,
                        help="Path to pickled dict with aminoacid2dna translations")
    parser.add_argument("--convert_2_amino", type=bool, default=False,
                        help="True if sequencing aminoacids")
    parser.add_argument("--convert_2_dna", type=bool, default=False,
                        help="True if data in RNA")
    parser.add_argument("--aligner_mode", type=str, default='local',
                        help="local or global.")
    args = parser.parse_args()

    aligner = Aligner([args.seq_A_path, args.seq_A_path], args.aminoacid_dict_path, args.aligner_mode, False, False)
    aligner.read_seqs()

    print(aligner.get_score())
    aligner.print_alignment()
