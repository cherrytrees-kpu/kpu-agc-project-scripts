"""
Generate a taxid map file for making a custom BLAST db
"""
from Bio import SeqIO
import pathlib
import argparse

def parse_args():
    parser = argparse.ArgumentParser(description='Generate a taxid map from .gb file data')
    parser.add_argument('f', 
        metavar='input', 
        action='store', 
        type=pathlib.Path, 
        help = 'name of input .gb file'
    )
    args = parser.parse_args()
    return args.f

def main():
    genbank_filepath = parse_args()
    num_seq = 0
    num_unidentified = 0
    #Open Genbank file
    output_path = genbank_filepath.with_name(genbank_filepath.stem + '_taxid_map.txt')
    f = open(output_path,'w')
    #Map db_xref to accesion
    for gb in SeqIO.parse(genbank_filepath,"gb"):
        try:
            annotations = gb.annotations['accessions'][0]
            for db_xref in gb.features[0].qualifiers['db_xref']:
                if db_xref.split(':')[0] == 'taxon': 
                    taxid = db_xref.split(':')[1]
            #taxid = gb.features[0].qualifiers['db_xref'][0].split(':')[1]
            f.write("{} {}\n".format(annotations, taxid))
            num_seq = num_seq + 1
        except:
            num_unidentified = num_unidentified + 1
            pass
    f.close()
    print(f'Total number of sequences: {str(num_seq)}')
    print(f'Total number of sequences with no species: {str(num_unidentified)}')

if __name__ == "__main__":
    main()