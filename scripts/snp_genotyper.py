from Bio import AlignIO, SeqIO
import argparse
import subprocess
import tempfile
import csv
import pathlib


def parse_args(): 
    """
    Arguments: 
    1) reference_path: path to the reference sequence fasta file
    2) ab1_path: path to all of the ab1 files
    """
    parser = argparse.ArgumentParser(description='Quick genotyper. Uses abi-trim to trim the ab1 files, and aligns against a reference sequence.')
    parser.add_argument(
        'reference_path',
        action='store',
        type=pathlib.Path, 
        help='Path to reference sequence' 
    )
    parser.add_argument(
        'ab1_path',
        action='store',
        type=pathlib.Path,
        help='Path to ab1 files'
    )
    args = parser.parse_args()
    return (args.reference_path, args.ab1_path)

def get_sequence_regions(alignment): 
    """ For each sequence in an alignment, determine its
    start and end alignment coordinates
    Parameters: 
    alignment - AlignIO object for the alignment
    Output: 
    list_regions - dictionary with start and end coordinates
    [accession id]:(start, end)
    """
    list_regions = dict()
    for accession in alignment: 
        #Variables
        f_start = False
        f_end = False
        start = 0
        end = 0
        #Data from SeqRecord object 
        sequence = accession.seq
        accession_id = (accession.id).split('.')[0]
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
        list_regions[accession_id] = (start, end)
    return list_regions

def generate_alignment(reference_seq, ab1_seqs): 
    alignment_input = dict()
    gene_name = reference_seq.id
    alignment_input[gene_name] = str(reference_seq.seq)
    for ab1_seq in ab1_seqs: 
        cow_id = ab1_seq.id.split('_')[1]
        alignment_input[cow_id] = str(ab1_seq.seq)
    #Generate temp fasta file
    fasta = tempfile.NamedTemporaryFile(delete=True)
    for seq_id in alignment_input: 
        fasta.write(f'>{seq_id}\n{alignment_input[seq_id]}\n'.encode())
    fasta.seek(0)
    args = [
        'mafft',
        '--auto',
        fasta.name,
    ]
    result = subprocess.run(args, capture_output=True)
    decoded = result.stdout.decode('utf-8')
    fasta.close()
    output_file = open(f'{gene_name}-aligned.fasta', 'w')
    output_file.write(decoded)
    output_file.close()

def main(reference_path, ab1_path): 
    """
    Workflow: 
    1) Open reference fasta file
    2) Open each of the ab1 files, and then import them into a list --> use abi-trim
    3) Using mafft, align all of the ab1 sequences to the reference fasta file. 
    4) Call the SNPs for each position of the alignment
    5) 
    """
    #Open reference sequence
    reference_seq = SeqIO.read(reference_path, 'fasta')
    #Open the ab1 files
    ab1_seqs = []
    for ab1_file_path in ab1_path.glob('*.ab1'): 
        sequence = SeqIO.read(ab1_file_path, 'abi-trim')
        ab1_seqs.append(sequence)
    #Create alignment
    generate_alignment(reference_seq, ab1_seqs)
    #Open Alignment
    alignment = AlignIO.read(f'{reference_seq.id}-aligned.fasta', 'fasta')
    #Define edges of the alignment
    list_regions = get_sequence_regions(alignment)
    del list_regions[reference_seq.id]
    start_pos = max([list_regions[x][0] for x in list_regions])
    end_pos = min([list_regions[x][1] for x in list_regions])
    #Datastructure for sequence data: 
    snp_header = []
    snp_header.append('Position')
    for seq in alignment: 
        snp_header.append(seq.id)    
    snp_data = []
    for position in range(start_pos, end_pos): 
        if len(set(alignment[:,position])) > 1:
            position_data = list(alignment[:,position])
            position_data.insert(0, position+1)
            snp_data.append(position_data)
    output_file = open('snp_calls.csv', 'w', newline='')
    csvwriter = csv.writer(output_file, dialect='excel')
    csvwriter.writerow(snp_header)
    csvwriter.writerows(snp_data)
    output_file.close()

if __name__ == '__main__': 
    reference_path, ab1_path = parse_args()
    main(reference_path, ab1_path)
