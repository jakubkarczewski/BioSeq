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
    parser.add_argument("--aligner_mode", type=str, default='global',
                        help="local or global.")
    parser.add_argument("--distance_matrix", type=str, default='distance_matrix.csv',
                        help="Distance matrix path.")
    args = parser.parse_args()

    aligner = Aligner([args.seq_A_path, args.seq_B_path], args.aminoacid_dict_path, args.aligner_mode,
                      args.distance_matrix, False, False)
    # read
    aligner.read_seqs()

    # rna 2 dna
    if args.convert_2_dna:
        aligner.rna_2_dna()

    # get alignments with aminoacid translation
    alignments = aligner.get_alignments(custom_matrix=True, analyse_amino=args.convert_2_amino)

    # get alignment with min score
    alignment_str = None
    min_score = None
    first_iter = True
    for alignment in alignments:
        if first_iter:
            min_score = alignment.score
            alignment_str = str(alignment)

        if min_score > alignment.score:
            min_score = alignment.score
            alignment_str = str(alignment)

    # print alignment info
    print('Alignment:\n', alignment_str)
    print('Similarity:', min_score)
    print('Edit distance:', 1 - min_score)
    print('Levenstein distance:', aligner.get_levenstein_distance(analyse_amino=args.convert_2_amino))
