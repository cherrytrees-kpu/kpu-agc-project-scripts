import argparse
import subprocess
import pathlib

def run_tracy(gene_name, sample_id, reference_fasta_path, ab1_file_path, peak_ratio, trim): 
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
        trim,
        '-p',
        peak_ratio,
        '-o',
        f'{gene_name}_{sample_id}',
        ab1_file_path,
    ]
    result = subprocess.run(args, capture_output=False)

def parse_args(): 
    parser = argparse.ArgumentParser(description='Batch execute Tracy')
    parser.add_argument('reference_file_path', 
        action='store', 
        type=pathlib.Path,
        help = 'Path to the reference sequence, fasta format'
    )
    parser.add_argument('ab1_path',
        action='store', 
        type=pathlib.Path, 
        help = 'Path to the directory where all of the ab1 files are stored'
    )
    parser.add_argument('-p',
        dest='peak_ratio',
        default='0.5',
        action='store',
        type=str,
        help='Peak ratio to pass to tracy. 0 to 1.'
    )
    parser.add_argument('-t',
        dest='trim',
        default='2',
        action='store',
        type=str,
        help='Trimming stringency to pass to Tracy. 1 to 9. 1 is most stringent.'
    )
    args = parser.parse_args()

    return (args.reference_file_path, args.ab1_path, args.peak_ratio, args.trim)

def main(reference_file_path, ab1_path, peak_ratio, trim): 
    gene_name = reference_file_path.stem.split('.')[0]
    for ab1_file_path in ab1_path.glob('*.ab1'): 
        cow_name = ab1_file_path.stem.split('_')[2]
        run_tracy(gene_name, cow_name, reference_file_path, ab1_file_path, peak_ratio, trim)

if __name__ == '__main__': 
    reference_file_path, ab1_path, peak_ratio, trim = parse_args()
    main(reference_file_path, ab1_path, peak_ratio, trim)

