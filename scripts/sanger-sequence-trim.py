import argparse
import pathlib
from Bio import SeqIO

def parse_args(): 
    parser = argparse.ArgumentParser(description='Process ab1 files and get Mott algorithm-trimmed sequences')
    parser.add_argument('sequences', 
        metavar='path_to_ab1_files', 
        action='store', 
        type=str,
        help = 'path to ab1 sequence files'
    )
    args = parser.parse_args()
    ab1_path = pathlib.Path(args.sequences)
    return ab1_path
def main(ab1_path):
    output_path = ab1_path.joinpath('output.fasta')
    output_file = open(output_path, 'w')

    ab1_iter = ab1_path.glob('*.ab1')
    for ab1_file_path in ab1_iter: 
        ab1_file = SeqIO.read(ab1_file_path, 'abi-trim')
        output_file.write(f'>{ab1_file.id}\n')
        output_file.write(f'{str(ab1_file.seq)}\n')
    output_file.close()

if __name__ == '__main__':
    ab1_path = parse_args()
    main(ab1_path)