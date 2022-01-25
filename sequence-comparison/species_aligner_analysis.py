from Bio import SeqIO, AlignIO
import argparse
import pathlib
import csv
import subprocess

class genbankHandler(): 
    def __init__(self, gb_path): 
        self.data = SeqIO.parse(gb_path, 'gb')
        self.path = pathlib.Path(gb_path)
        self.fasta_paths = []
        self.alignment_paths = []
        self.gb = dict()
        self.species_keys = self.gb.keys()
        self._split_gb_by_species()
    def _split_gb_by_species(self):
        for gb_data in self.data: 
            organism = gb_data.annotations['organism']
            species_name = organism.split(' ')[1]
            #If there's no clear species designation, ignore
            if species_name != 'sp.':
                if organism not in self.gb.keys(): 
                    self.gb[organism] = [gb_data]
                else: 
                    self.gb[organism].append(gb_data)
    def num_species(self): 
        return len(self.gb.keys()) 
    def output_species_gb(self): 
        for organism in self.gb.keys(): 
            species_path = self.path.with_name(f'{organism}_data.gb')
            SeqIO.write(self.gb[organism], species_path, 'gb')
    def output_species_fasta(self): 
        for organism in self.gb.keys(): 
            species_path = self.path.with_name(f'{organism}.fasta')
            SeqIO.write(self.gb[organism], species_path, 'fasta')
            self.fasta_paths.append(pathlib.Path(species_path))
    def output_metadata(self): 
        #Loop through each organism
        for organism in self.gb.keys(): 
            species_path = self.path.with_name(f'{organism}_metadata.csv')
            all_gb_data = []
            #Loop through each Genbank file
            for gb in self.gb[organism]: 
                gb_metadata = [
                    gb.id, 
                    len(gb.seq), 
                    gb.description,
                    gb.annotations['source'], 
                ]
                #Identify source feature, then store data
                for feature in gb.features: 
                    if feature.type == 'source': 
                        if 'isolate' in feature.qualifiers: 
                            gb_metadata.append(feature.qualifiers['isolate'])
                        else: 
                            gb_metadata.append('')
                        if 'host' in feature.qualifiers: 
                            gb_metadata.append(feature.qualifiers['host'])
                        else: 
                            gb_metadata.append('')
                        if 'country' in feature.qualifiers: 
                            gb_metadata.append(feature.qualifiers['country'])
                        break
                all_gb_data.append(gb_metadata)
            metadata_output = open(species_path, 'w', newline='')
            csv_writer = csv.writer(metadata_output)
            csv_writer.writerow(
                (
                    'accession_id',
                    'seq_len',
                    'description',
                    'organism',
                    'isolate',
                    'host',
                    'country'
                )
            )
            csv_writer.writerows(all_gb_data)
            metadata_output.close()
    def generate_species_alignment(self): 
        for fasta_path in self.fasta_paths: 
            aligned_path = fasta_path.with_name(fasta_path.stem + '_aligned.fasta')
            args =[
                'mafft',
                '--auto',
                str(fasta_path), 
            ]
            result = subprocess.run(args, capture_output=True)
            decoded = result.stdout.decode('utf-8')
            output_file = open(aligned_path, 'w')
            output_file.write(decoded)
            output_file.close()
            self.alignment_paths.append(pathlib.Path(aligned_path))
    def generate_consensus_sequence(self): 
        def _get_consensus(): 
            pass
        for alignment_path in self.alignment_paths: 
            alignment = AlignIO.parse(alignment_path, 'fasta')


def parse_args(): 
    parser = argparse.ArgumentParser('Analyze and align multiple sequences for genus-level analysis')
    parser.add_argument(
        'genus_gb', 
        action='store',
        type=pathlib.Path,
        help='path to genbank file containing all data'
    )
    args = parser.parse_args()
    return args.genus_gb

def main(): 
    """
    Algorithm: 
    1) Open genus genbank file --> this is generated using Entrez query on NCBI Nucleotide
    2) Split gb into several species genbank files
    3) Align .fasta files for each species

    """
    genus_gb = parse_args()
    genus_data = genbankHandler(genus_gb)
    genus_data.output_species_gb()
    genus_data.output_species_fasta()
    genus_data.output_metadata()
    genus_data.generate_species_alignment()

if __name__ == '__main__': 
    main()