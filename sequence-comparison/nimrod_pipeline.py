from Bio import SeqIO
import argparse
import pathlib
import pandas
import subprocess

def split_gb_by_species(gb_path):
    """
    Split a .gb into multiple gb entries by species
    """
    #ALGORITHM
    #1) Open .gb file
    data = SeqIO.parse(gb_path, 'gb')
    #Store all data into list
    list_gb = []
    for gb in data: 
        list_gb.append(gb)
    species_data = dict()
    #Loop will store all species data into a dictionary, where key provides list of gb entries
    for gb in list_gb: 
        #Process name
        organism = gb.annotations['organism']
        #e.g. Meloidogyne spartelensis --> 2nd element is species name
        species_name = organism.split(' ')[1]
        #If there's no clear species designation, ignore
        if species_name != 'sp.':
            if organism not in species_data.keys(): 
                species_data[organism] = [gb]
            else: 
                species_data[organism].append(gb)
    return species_data
def get_species_seq_count(species_data): 
    """
    Count number of species within the genus, and list the number of sequences for that species. 
    """
    species_seq_count = []
    for species in species_data.keys(): 
        species_seq_count.append((species, len(species_data[species])))
    species_seq_count.sort(reverse=True, key=lambda x:x[1])
    return species_seq_count
def output_species_gb(species_data, gene_name, output_path):
    """
    Write the Genbank file for a specific species
    """
    for species in species_data.keys():
        species_path_name = species.replace(' ', '_')
        species_path = output_path.joinpath(f'{gene_name}_{species_path_name}_data.gb')
        SeqIO.write(species_data[species], species_path, 'gb')
    #Output combined
    all_species_data = []
    genus_name = next(iter(species_data.values()))[0].annotations['organism'].split(' ')[0]
    for species in species_data.keys(): 
        all_species_data.extend(species_data[species])
    all_species_path = output_path.joinpath(f'{gene_name}_{genus_name}_data.gb')
    SeqIO.write(all_species_data, all_species_path, 'gb')
def output_metadata(species_seq_count, gene_name, output_path):
    """
    Provide information on the final dataset. 
    """
    genus_name = species_seq_count[0][0].split(' ')[0]
    metadata_file = open(output_path.joinpath(f'{gene_name}-{genus_name}-metadata.txt'), 'w')
    metadata_file.write(f'Total number of species: {str(len(species_seq_count))}\n')
    for species in species_seq_count: 
        metadata_file.write(f'{species[0]}\t{str(species[1])}\n')
    metadata_file.close()
def gb_extract(gb_path):
    def _retrieve_metadata(gb_entry):
        """Parse GenBank entry to obtain relevant metadata information
        Parameters: 
        gb_entry - SeqRecord object
        Returns: 
        metadata - tuple of metadata
        """ 
        def get_source_index(features): 
            """Get the index of the 'source' entry for GenBank entry
            Parameters: 
            features - feature table of a GenBank entry
            Returns: 
            index - int - position of 'source' entry in feature table
            """
            index = None
            for feature in features: 
                if feature.type == "source":
                    index = features.index(feature)
            return index
        def get_source_qualifiers(source_qual): 
            """Get qualifiers of interest for source
            Parameters: 
            source_qual - qualifiers in source
            Returns: 
            qual - list - contains qualifiers of interest
            """
            quals = []
            interest = ('isolate', 'host', 'country')
            for qual in interest: 
                if qual in source_qual:
                    if len(source_qual[qual]) == 1: 
                        quals.append(source_qual[qual][0])
                    else: 
                        quals.append(','.join(source_qual[qual]))
                else: 
                    quals.append(None)
            return quals
        #Get qualifiers from 'source' feature
        source_index = get_source_index(gb_entry.features)
        source_qualifiers = get_source_qualifiers(gb_entry.features[source_index].qualifiers)
        #Order: accession id, sequence length, description, organism, isolate, host, country
        metadata = (
            gb_entry.id, 
            len(gb_entry.seq), 
            gb_entry.description, 
            gb_entry.annotations['source'],
            source_qualifiers[0],
            source_qualifiers[1],
            source_qualifiers[2],
        )
        return metadata
    def _output_fasta(list_fasta, path):
        """Create .fasta files
        Parameters: 
        list_fasta - list - list containing fasta strings
        path - pathLib.path - path without suffix
        Returns: 
        None
        """
        output_file = open(path, 'w')
        for fasta in list_fasta: 
            output_file.write(fasta)
        output_file.close()
    data = SeqIO.parse(gb_path, 'gb')
    list_fasta = []
    metadata = {
            'acc':[],
            'seq_len':[],
            'desc':[],
            'organism':[],
            'isolate':[],
            'host':[],
            'country':[],
        }
    for entry in data: 
        list_fasta.append(entry.format('fasta'))
        meta = _retrieve_metadata(entry)
        #Add to metadata
        metadata['acc'].append(meta[0])
        metadata['seq_len'].append(meta[1])
        metadata['desc'].append(meta[2])
        metadata['organism'].append(meta[3])
        metadata['isolate'].append(meta[4])
        metadata['host'].append(meta[5])
        metadata['country'].append(meta[6])
        meta_df = pandas.DataFrame(data=metadata)
        #Output
        _output_fasta(list_fasta, gb_path.with_name(gb_path.stem + '.fasta'))
        meta_df.to_csv(gb_path.with_name(gb_path.stem + '_metadata.csv'), index=False)
def execute_mafft(fasta_path): 
    """
    Execute MAFFT using subprocess
    Manual: 
    https://mafft.cbrc.jp/alignment/software/manual/manual.html
    """
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
def parse_args(): 
    parser = argparse.ArgumentParser(description='Generate .fasta files for each unique species. Ignores genus-only designations.')
    parser.add_argument('d', metavar='input', type=str, help = 'path to a GenBank file')
    parser.add_argument('gene',metavar='gene_Name', type=str, help='name of target site')
    args = parser.parse_args()
    return (pathlib.Path(args.d), args.gene)
def main(gb_path, gene_name):
    #Algorithm: 
    #1) Open genus genbank file
    #2) Split genus genbank into multiple genbank files
    #3) Output these files, and output metadata
    #4) For each species genbank file, extract appropriate data
    #5) For each species .fasta file that has more than one sequence, 
    #   align using MAFFT 
    output_path = gb_path.parent.joinpath('output')
    output_path.mkdir(exist_ok=True)
    species_data = split_gb_by_species(gb_path)
    species_seq_count = get_species_seq_count(species_data)
    output_species_gb(species_data, gene_name, output_path)
    output_metadata(species_seq_count, gene_name, output_path)
    species_gb_iter = output_path.glob('*.gb')
    for species_gb_path in species_gb_iter: 
        gb_extract(species_gb_path)
    species_fasta_iter = output_path.glob('*.fasta')
    for species_fasta_path in species_fasta_iter: 
        execute_mafft(species_fasta_path)

if __name__ == '__main__': 
    gb_path, gene_name = parse_args()
    main(gb_path, gene_name)