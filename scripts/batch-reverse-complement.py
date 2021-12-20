from Bio import Seq, SeqIO
import argparse
import pathlib

def parse_args(): 
    parser = argparse.ArgumentParser(description="Script to quickly generate reverse complement fastas")
    parser.add_argument(
        'seq_path',
        type=pathlib.Path,
        action='store',
        help='Path to sequence files'
    )
    args = parser.parse_args()
    return args.seq_path
def main(seq_path):
    if seq_path.is_dir() is True: 
        #Loop over all fasta files in the directory
        input_paths = []
        for seq_file_path in seq_path.glob('*.fasta'): 
            input_paths.append(seq_file_path)
        for seq_file_path in input_paths:
            seq_file = SeqIO.read(seq_file_path, 'fasta')
            gene = seq_file.id
            rev_comp_seq = seq_file.seq.reverse_complement()
            output_file_path = seq_file_path.with_name(f'{gene}_REV.fasta')
            output_file = open(output_file_path, 'w')
            output_file.write(
                f'>{gene}\n{str(rev_comp_seq)}'
            )
            output_file.close()
    else: 
        seq_file = SeqIO.read(seq_path, 'fasta')
        gene = seq_file.id
        rev_comp_seq = seq_file.seq.reverse_complement()
        output_file_path = seq_path.with_name(f'{gene}_REV.fasta')
        output_file = open(output_file_path, 'w')
        output_file.write(
            f'>{gene}\n{str(rev_comp_seq)}'
        )
        output_file.close()

if __name__ == '__main__': 
    seq_path = parse_args()
    main(seq_path)