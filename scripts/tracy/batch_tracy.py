from Bio import SeqIO
import argparse
import subprocess
import pathlib

def find_trim_values(ab1_file_path):
    untrimmed = SeqIO.read(ab1_file_path, 'abi')
    trimmed = SeqIO.read(ab1_file_path, 'abi-trim')
    left_trim = untrimmed.seq.find(trimmed.seq[0:5])-1
    right_trim = len(untrimmed.seq) - len(trimmed) - left_trim
    return left_trim, right_trim

def run_tracy(gene_name, sample_id, reference_fasta_path, ab1_file_path, peak_ratio, left_trim, right_trim): 
    """
    Run tracy using the reference fasta files
    """
    args=[
        'tracy',
        'decompose',
        '-r',
        reference_fasta_path,
        '-v',
        '-t',
        '0',
        '-p',
        peak_ratio,
        '-o',
        f'{gene_name}_{sample_id}',
        '-q',
        str(left_trim),
        '-u',
        str(right_trim), 
        ab1_file_path,
    ]
    result = subprocess.run(args, capture_output=False)

def parse_args(): 
    parser = argparse.ArgumentParser(description='Batch execute Tracy')
    parser.add_argument('reference_file_path', 
        action='store', 
        type=pathlib.Path,
        help = 'Path to the reference sequence/genome'
    )
    parser.add_argument('ab1_path',
        action='store', 
        type=pathlib.Path, 
        help = 'Path to the directory where all of the ab1 files are stored'
    )
    parser.add_argument('-p',
        dest='peak_ratio',
        default='0.8',
        action='store',
        type=str,
        help='Peak ratio to pass to tracy. 0 to 1.'
    )

    args = parser.parse_args()

    return (args.reference_file_path, args.ab1_path, args.peak_ratio)

def main(): 
    
    reference_file_path, ab1_path, peak_ratio = parse_args()
    
    for ab1_file_path in ab1_path.glob('*.ab1'):
        gene_name = ab1_file_path.stem.split('_')[1] 
        cow_name = ab1_file_path.stem.split('_')[2]
        left_trim, right_trim = find_trim_values(ab1_file_path)
        run_tracy(gene_name, cow_name, reference_file_path, ab1_file_path, peak_ratio, left_trim, right_trim)

if __name__ == '__main__': 
    main()