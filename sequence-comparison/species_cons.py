"""
species_cons.py - pipeline to generate consensus sequences from species alignments, 
then align these consensus sequences to each other for creating sequence phylogenetic
trees
"""

from Bio import AlignIO, Seq
import csv
import argparse
import subprocess
import pathlib

def get_seq_position_info(alignment): 
    """
    Get information for each position of the alignment.
    Output is in list format, with each index corresponding to the 0-based position
    of the alignment. 
    Each index contains dictionary with following information: 
    pos - position number
    accessions - the accessions represented at that position
    bases - base calls in corresponding order of accessions
    seqs_rep - number of accessions represented at that position
    p_seq_rep - percentage of accessions represented at that position
    """
    def get_sequence_regions(alignment): 
        """
        Identify the regions that each accession spans on the alignment.
        Algorithm: 
        1) From beginning of sequence, keep going until position is not '-',
           then record the position as the start of the region
        2) From the end of the sequence, go in reverse until position is not '-',
           then record the position as the end of the region
        3) Return tuple --> (accession, start, end)
        """
        def get_sequence_region(accession): 
            #Variables
            f_start = False
            f_end = False
            start = 0
            end = 0
            #Data from SeqRecord object 
            sequence = accession.seq
            acc = (accession.id).split('.')[0]
            #From beginning of the sequence, keep going until position is not '-'
            i = 0
            while (f_start is False): 
                if sequence[i] != '-': 
                    start = i
                    f_start = True
                i = i + 1
            #From ending of the sequence, keep going until position is not '-'
            i = -1
            while (f_end is False): 
                if sequence[i] != '-':
                    end = i%len(sequence) #Translate negative index to positive index
                    f_end = True
                i = i - 1
            return ((acc, start, end))
        sequence_regions = []
        for accession in alignment: 
            sequence_regions.append(get_sequence_region(accession))
        return sequence_regions
    def get_contributing_seqs(position, list_seq_region):
        """
        Given a position and the list of sequence regions for each accession,
        identify which accessions are represented
        """
        list_accession = []
        list_indices = []
        for sequence in list_seq_region:
            if (position >= sequence[1]) and (position <= sequence[2]): 
                list_accession.append(sequence[0])
                list_indices.append(list_seq_region.index(sequence))
        return list_accession, list_indices
    list_seq_position_info = []
    #Get sequence regions
    list_seq_regions = get_sequence_regions(alignment)
    num_accessions = len(alignment)
    for position in range(alignment.get_alignment_length()):
        list_accessions, list_indices = get_contributing_seqs(position, list_seq_regions)
        list_bases = []
        #Append the basecalls of every contributing sequence at specified position
        for index in list_indices: 
            list_bases.append(alignment[index].seq[position])
        position_info = {
            'pos':position,
            'accessions':tuple(list_accessions),
            'bases':tuple(list_bases),
            'seqs_rep':len(list_accessions),
            'p_seq_rep':len(list_accessions)/num_accessions,
        }
        list_seq_position_info.append(position_info)
    return list_seq_position_info
def get_alignment_summary(seq_info): 
    """
    Determine the consensus sequence of an alignment, and create position matrix
    Definition of consensus: most common base represented at that position. 
    """
    consensus_sequence = []
    position_matrix = []
    for position in seq_info: 
        #Ignore any ambiguous basecalls - accept A, T, C, G, and 'gap'
        base_counts = {
            'a':position['bases'].count('a')+position['bases'].count('A'),
            't':position['bases'].count('t')+position['bases'].count('T'),
            'c':position['bases'].count('c')+position['bases'].count('C'),
            'g':position['bases'].count('g')+position['bases'].count('G'),
            '-':position['bases'].count('-'),
        }
        #print(base_counts)
        max_basecalls = [key for key, count in base_counts.items() if count == max(base_counts.values())]
        if len(max_basecalls) == 1: 
            consensus_sequence.append(max_basecalls[0])
        else: 
            consensus_sequence.append('n')
        #Assembling position_matrix
        position_matrix.append(base_counts)
    return (''.join(consensus_sequence), position_matrix)
def output_csv(data, output_path, name): 
    output_file_path = output_path.joinpath(name)
    output_file = open(output_file_path, 'w', newline='')
    csv_dict_writer = csv.DictWriter(output_file, data[0].keys())
    csv_dict_writer.writeheader()
    csv_dict_writer.writerows(data)
    output_file.close()
def create_fasta(cons_seqs, output_path, name): 
    output_file_path = output_path.joinpath(name)
    output_file = open(output_file_path, 'w')
    for key in cons_seqs: 
        output_file.write(f'>{key}\n')
        output_file.write(cons_seqs[key])
        output_file.write('\n')
    output_file.close()
def parse_args():
    parser = argparse.ArgumentParser(description='species_cons.py - pipeline to generate species consensus sequences')
    parser.add_argument('species_alignment_path', 
        metavar='species_alignment', 
        action='store', 
        type=pathlib.Path,
        help = 'path to species alignments directory'
    )
    parser.add_argument('-b, --batch',
        dest='flag_batch',
        action='store_true', 
        help='Use flag if you want to pass a directory'
    )
    args = parser.parse_args()

    return args.species_alignment_path, args.flag_batch
def main():

    species_alignment_path, flag_batch = parse_args()

    if flag_batch is False:
        output_path = species_alignment_path.parent 
        alignment = AlignIO.read(species_alignment_path, 'fasta')
        alignment_id = species_alignment_path.stem
        seq_info = get_seq_position_info(alignment)
        cons_seq, position_matrix = get_alignment_summary(seq_info)
        #UNGAP
        ungapped_cons_seq = str(Seq.Seq(cons_seq).ungap())
        #Output position_matrix
        output_csv(position_matrix, output_path, alignment_id+'_pos-matrix.csv')
        #Output seq_info
        output_csv(seq_info, output_path, alignment_id+'_seq-info.csv')
        #Output consensus sequence
        fasta_path = output_path.joinpath(f'{alignment_id}_consensus.fasta')
        fasta_file = open(fasta_path, 'w')
        fasta_file.write(f'>{alignment_id}\n')
        fasta_file.write(ungapped_cons_seq)
        fasta_file.close()
    else: 
        return   

    #output_path = species_alignment_path.joinpath('results')
    #output_path.mkdir(exist_ok=True)
    #cons_seqs = dict()
    #for alignment_path in species_alignment_path.glob('*.fasta'):
    #    alignment = AlignIO.read(alignment_path, 'fasta')
    #    alignment_id = alignment_path.name.split('_data')[0]
    #    seq_info = get_seq_position_info(alignment)
    #    cons_seq, position_matrix = get_alignment_summary(seq_info)
        #Output position_matrix
    #    output_csv(position_matrix, output_path, alignment_id+'_pos-matrix.csv')
        #Output seq_info
    #    output_csv(seq_info, output_path, alignment_id+'_seq-info.csv')
    #    cons_seqs[alignment_id] = cons_seq
    #create_fasta(cons_seqs, output_path, 'all_sequences.fasta')

if __name__ == '__main__':

    main()