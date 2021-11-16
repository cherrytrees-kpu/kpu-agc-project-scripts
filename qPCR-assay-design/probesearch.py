from Bio import SeqIO, Seq
import pathlib
import argparse
import tempfile
import io
import subprocess
import pandas
import csv

MIN_PRIMER_LEN = 17
MAX_PRIMER_LEN = 22
LEN_DB = 16194
BLASTDB = 'SSU_nematode.fasta'

def check_primer(seq): 
    def count_c_g(seq): 
        return seq.count('c') + seq.count('g') + seq.count('C') + seq.count('G')
    def chk_run(seq): 
        """
        Ensure that there are no runs of four or more identical nucleotides in the probe
        """
        current = seq[0]
        identical_len = 1
        detected = False
        #Go through each nucleotide in the primer
        for nucl in range(1, len(seq)): 
            if seq[nucl] == current: 
                identical_len = identical_len + 1
                if identical_len > 3: 
                    detected = True
                    break
            else: 
                current = seq[nucl]
                identical_len = 1
        return detected
    def chk_last_5(seq): 
        return count_c_g(seq[-5:]) > 2
    #Check the three conditions
    percent_c_g = count_c_g(seq)/len(seq)
    run_flag = chk_run(seq)
    last5_flag = chk_last_5(seq)
    #print(f"{str(percent_c_g)}{str(run_flag)} {str(last5_flag)}")
    #Return true if all of the conditions are passed, otherwise return false
    if(
        0.3 <= percent_c_g <= 0.8 
        and run_flag is False
        and last5_flag is False
        and seq[0] != 'G'
        and seq[0] != 'g'
    ):
        #print('Passed!')
        return True
    else: 
        return False

def get_probes(target_sequence, target_start): 
    """
    Function finds all viable probes. 
    Probe criteria: 
    1) 5' end of the probe is not a G
    2) More C's than G's
    3) No more than 4 nucleotides in a row
    4) 30 - 80% CG content
    """
    list_probes = []
    #For each root position in the target sequence
    for i in range(len(target_sequence) - MIN_PRIMER_LEN + 1):
        #Iterate over all probe lengths
        for probe_len in range(MIN_PRIMER_LEN, MAX_PRIMER_LEN + 1): 
            probe_seq = target_sequence[i:i+probe_len]
            if check_primer(probe_seq) is True: 
                list_probes.append(
                    {
                        'root_pos':target_start + i,
                        'len':probe_len,
                        'seq':str(probe_seq),
                    }
                )
    return list_probes

def generate_blast_results_db(list_probes): 
    blast_results = dict()
    for probe in list_probes:
        probe_id = f"{str(probe['root_pos'])}-{probe['len']}" 
        blast_results[probe_id] = blast(probe['seq'], BLASTDB, LEN_DB)
    return blast_results

def blast(seq, blastdb, num_aligned): 
    fasta = tempfile.NamedTemporaryFile(delete=True)
    fasta.write(f">probe\n{str(seq)}".encode())
    fasta.seek(0)
    args = [
        "blastn",
        "-task",
        "blastn-short",
        "-db",
        blastdb,
        "-num_alignments",
        str(num_aligned),
        "-outfmt",
        "10 qacc sacc ssciname pident qlen length mismatch gapopen qstart qend sstart send evalue bitscore",
        "-query",
        fasta.name,
    ]
    #Capture output
    result = subprocess.run(args, capture_output=True)
    decoded = result.stdout.decode('utf-8')
    print(decoded)
    output = io.StringIO(decoded)
    #Output formatting into dataframe
    headers=[
        'qacc',
        'sacc',
        'ssciname',
        'pident',
        'qlen',
        'length',
        'mismatch', 
        'gapopen', 
        'qstart', 
        'qend', 
        'sstart', 
        'send', 
        'evalue', 
        'bitscore',
    ]
    data = pandas.read_csv(output, sep=',', header=None, names=headers)
    fasta.close()
    return data

def calculate_sensitivity(blast_results, target_accessions):
    """
    Function calculates sensitivity of the oligo (binding to target sequences). 
    Binding is defined as 100% query coverage and 100% percent identity of the oligo to the target sequence.
    Sensitivity is calculated by the following formula: 
        TP / (TP + FN)
        where: 
        TP = succesfully amplified accessions
        FN = target accessions that were not amplified
        TP + FN = total number of target accessions
    """
    #Take only accessions where there was a perfect match --> full query coverage, 100% identity
    perfect_match_results = blast_results.loc[(blast_results['qlen']==blast_results['length'])&(blast_results['pident']==100.0)]
    #Retrieve only the accessions list
    amplified_accessions = set(perfect_match_results.loc[:,'sacc'])
    target_match = 0
    #Count number of target accessions in the amplified accession list 
    for accession in amplified_accessions: 
            if accession in target_accessions: 
                target_match = target_match + 1
    #Calculate sensitivity and return
    sensitivity = target_match/len(target_accessions)
    return sensitivity

def calculate_specificity(blast_results, target_accessions, blastdb_len): 
    """
    Function calculates specificity of the oligo (binding to non-target sequences). 
    Binding is defined as simply appearing in the BLAST results --> this will be a overestimation of the specificity,
    making it the 'worst case' scenario. 
    Specificity is calculated by the following formula: 
        TN / (TN + FP)
        where:
        TN = non-target accessions that were not amplified
        TN = (total non-target) - (amplified non-target)
        FP = amplified non-target
        resulting in the following formula: 
        ((Total non-target) - (amplified non-target)) / total non-target
    """
    blast_match_accessions = set(blast_results.loc[:,'sacc'])
    #Remove every target_accession from the blast_match_accessions list
    for accession in target_accessions:
        if accession in blast_match_accessions:  
            blast_match_accessions.remove(accession)
    #Calculate specificity
    #Total non-target = all_blast_sequences - target_accessions
    specificity = (blastdb_len - len(blast_match_accessions))/(blastdb_len - len(target_accessions))
    return specificity

def parse_args(): 
    parser = argparse.ArgumentParser(description='probesearch.py - identify viable probes in an alignment for given target sequences')
    parser.add_argument('target_seq', 
        metavar='genus_alignments', 
        action='store', 
        type=str,
        help = 'path to species alignments directory'
    )
    parser.add_argument('target_accessions', 
        metavar='target_accessions', 
        action='store', 
        type=str,
        help = 'path to target_accessions'
    )
    parser.add_argument('target_start',
        metavar='target_start', 
        action='store',
        type=int,
        help='start coordinate of target region'
    )
    parser.add_argument('target_end',
        metavar='target_end', 
        action='store',
        type=int,
        help='end coordinate of target region'
    )
    args = parser.parse_args()
    target_seq_path = pathlib.Path(args.target_seq)
    target_accessions_path = pathlib.Path(args.target_accessions)
    return (target_seq_path, target_accessions_path, args.target_start, args.target_end)
def main(target_seq_path, target_accessions, target_start, target_end):
    headers=[
        'qacc',
        'sacc',
        'ssciname',
        'pident',
        'qlen',
        'length',
        'mismatch', 
        'gapopen', 
        'qstart', 
        'qend', 
        'sstart', 
        'send', 
        'evalue', 
        'bitscore',
    ]
    #Input data - alignment, target indices
    target_seq_file = SeqIO.read(target_seq_path, 'fasta')
    target_sequence = target_seq_file.seq[target_start:target_end]
    print(len(target_sequence))
    #Generate probes
    list_probes = get_probes(str(target_sequence), target_start)
    blast_results = generate_blast_results_db(list_probes)

    #Calculate sensitivity and specificity for each probe
    for probe in list_probes: 
        probe_id = f"{str(probe['root_pos'])}-{probe['len']}" 
        probe_blast_results = blast_results[probe_id]
        probe['sensitivity'] = calculate_sensitivity(probe_blast_results, target_accessions)
        probe['specificity'] = calculate_specificity(probe_blast_results, target_accessions, LEN_DB)
        probe['score'] = probe['sensitivity'] + probe['specificity']
    
    #Sort
    list_probes.sort(key=lambda x: x['score'], reverse=True)

    #Outputs:
    #1) blast_csv for each blast result for each primer
    #2) file listing all of the probes
    #Blast results
    blast_folder_path = target_seq_path.with_name('probe_blast')
    blast_folder_path.mkdir()
    for blast_result_key in blast_results: 
        blast_output_path = blast_folder_path.joinpath(f"{blast_result_key}_probe.csv")
        blast_output_file = open(blast_output_path, 'w')
        blast_results[blast_result_key].to_csv(blast_output_file)
        blast_output_file.close()
    #Probes
    probe_data = []
    for probe in list_probes: 
        probe_data.append(
            (
                probe['root_pos'],
                probe['len'], 
                probe['seq'],
                probe['sensitivity'],
                probe['specificity'],
                probe['score'],
            )
        )
    csv_output_path = target_seq_path.with_name('probe_candidates.csv')
    csvfile = open(csv_output_path, 'w', newline='')
    csvwriter = csv.writer(csvfile)
    csvwriter.writerow(
        (
            'probe_root',
            'probe_len',
            'probe_seq',
            'sens',
            'spec',
            'score',
        )
    )
    csvwriter.writerows(probe_data)
    csvfile.close()
if __name__ == '__main__': 
    target_seq_path, target_accessions_path, target_start, target_end = parse_args()
    target_accessions = []
    input_file = open(target_accessions_path, 'r')
    for line in input_file: 
        target_accessions.append(line.strip('\n'))
    input_file.close()
    main(target_seq_path, target_accessions, target_start, target_end)