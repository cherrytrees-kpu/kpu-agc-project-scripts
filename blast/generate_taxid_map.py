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
        type=str, 
        help = 'name of input .gb file'
    )
    args = parser.parse_args()
    return pathlib.Path(args.f)

def main(genbank_filepath):
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
        except:
            pass
    f.close()

if __name__ == "__main__":
    genbank_path = parse_args()
    main(genbank_path)