import argparse
import pathlib
import csv
from Bio import SeqIO

def parse_args(): 
    parser = argparse.ArgumentParser(description='Process ab1 files and get Mott algorithm-trimmed sequences')
    parser.add_argument(
        'ab1_path', 
        metavar='path_to_ab1_files', 
        action='store', 
        type=pathlib.Path,
        help = 'path to ab1 sequence files'
    )
    parser.add_argument(
        '-s', '--sep',
        dest='sep_flag',
        action='store_true',
        help='Flag to create separate files for each sequence'
    )
    args = parser.parse_args()
    return (args.ab1_path, args.sep_flag)

def main():
    ab1_path, sep_flag = parse_args()
    metadata = [] 
    if sep_flag is True:
        ab1_iter = ab1_path.glob('*.ab1')
        for ab1_file_path in ab1_iter:
            #Open
            ab1_file = SeqIO.read(ab1_file_path, 'abi')
            ab1_file_trim = SeqIO.read(ab1_file_path, 'abi-trim')
            #Check trims
            left_trim = ab1_file.seq.find(ab1_file_trim.seq[0:5])-1
            right_trim = len(ab1_file.seq) - len(ab1_file_trim) - left_trim
            metadata.append((ab1_file_trim.id, left_trim, right_trim))
            #Output
            output_path = ab1_file_path.with_name(f'{ab1_file_trim.id}_trimmed.fasta')
            SeqIO.write(ab1_file_trim, output_path, 'fasta')
    else: 
        output_path = ab1_path.joinpath('output.fasta')
        output_file = open(output_path, 'w')
        ab1_iter = ab1_path.glob('*.ab1')
        for ab1_file_path in ab1_iter:
            #Open
            ab1_file = SeqIO.read(ab1_file_path, 'abi') 
            ab1_file_trim = SeqIO.read(ab1_file_path, 'abi-trim')
            #Check trims
            left_trim = ab1_file.seq.find(ab1_file_trim.seq[0:5])-1
            right_trim = len(ab1_file.seq) - len(ab1_file_trim) - left_trim
            metadata.append((ab1_file_trim.id, left_trim, right_trim))
            #Output
            output_file.write(f'>{ab1_file_trim.id}\n')
            output_file.write(f'{str(ab1_file_trim.seq)}\n')
        output_file.close()
    metadata_path = ab1_path.joinpath('metadata.csv')
    metadata_file = open(metadata_path, 'w', newline='')
    csvwriter = csv.writer(metadata_file)
    csvwriter.writerow(('ID','left_trim','right_trim'))
    csvwriter.writerows(metadata)
    metadata_file.close()

if __name__ == '__main__':
    main()