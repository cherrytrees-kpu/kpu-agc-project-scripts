from Bio import SeqIO
import pathlib
import argparse

def parse_args(): 
    parser = argparse.ArgumentParser(description='Extract sequences from GB as fasta files')
    parser.add_argument(
        'd', 
        metavar='input', 
        type=pathlib.Path, 
        help = 'path to Genbank files')
    args = parser.parse_args()
    return args.d

def main(genbank_path):
    for genbank_file_path in genbank_path.glob('*.gb'): 
        #Read the Genbank file
        genbank = SeqIO.read(genbank_file_path, 'gb')
        gene_name = genbank_file_path.stem.split('_')[0]
        #Create fasta file, then write
        fasta_output_path = genbank_file_path.with_name(f'{gene_name}.fasta')
        fasta_file = open(fasta_output_path, 'w')
        fasta_file.write(f'>{gene_name}'+'\n')
        fasta_file.write(str(genbank.seq))
        fasta_file.close()

if __name__ == "__main__":
    genbank_path = parse_args()
    main(genbank_path)
