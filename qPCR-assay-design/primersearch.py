from Bio import SeqIO, Seq
import primer3
from operator import attrgetter
import pathlib
import argparse
import tempfile
import io
import subprocess
import pandas
import csv

MIN_PRIMER_LEN = 17
MAX_PRIMER_LEN = 22
#MAX_TM_DIFF = 5
LEN_DB = 16194
BLASTDB = 'SSU_nematode.fasta'

def parse_args(): 
    parser = argparse.ArgumentParser(description='Search for primers')
    parser.add_argument('target_seq_path', 
        metavar='target_seq', 
        action='store', 
        type=str, 
        help = 'path to target seq fasta'
    )
    parser.add_argument('list_target_accessions',
        metavar='target_accessions', 
        action='store', 
        type=str, 
        help = 'path to target accessions'
    )
    parser.add_argument('pb_start',
        metavar='pb_start', 
        action='store', 
        type=int, 
        help = 'start coordinate of probe'
    )
    parser.add_argument('pb_end',
        metavar='pb_end', 
        action='store', 
        type=int, 
        help = 'end coordinate of probe'
    )
    parser.add_argument('-d',
        dest='tm_diff',
        metavar='maximum_tm_diff',
        action='store',
        type=float,
        default=5.0,
        help='maximum temperature difference between forward and reverse'
    )
    args = parser.parse_args()
    target_seq_file_path = pathlib.Path(args.target_seq_path)
    target_accessions_path = pathlib.Path(args.list_target_accessions)
    pb_start = args.pb_start
    pb_end = args.pb_end
    return (target_seq_file_path, target_accessions_path, pb_start, pb_end, args.tm_diff)

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
    ):
        #print('Passed!')
        return True
    else: 
        return False

def find_fw_primers(target_seq, pb_start):
    """
    Function finds all viable FW primers.
    FW primer criteria: 
    1) 3'-end of the FW needs to be within 50 bp of the 5'-end of the probe.
    2) Between min and max length 
    3) % GC is 30% to 80%
    4) Last five nucleotides at the 3' end contain no more than two G + C residues
    5) No more than 4 consecutive nucleotides within the primer 

    Input data: 
    1) target_seq - str - target sequence - ATCTGATCATGATCATGACTAGTCATGGC
    2) pb_start - int - start index of 5'-end of the probe - 607

    Output data: 
    [
        {root_pos:0, len:17, seq:''}
    ]

    Algorithm: 
    1st root position is the nucleotide right before the 5'-end of the probe. 
    1st_root_pos_index = pb_start - 1
    However, for list-slicing purposes, take the pb_start index. 
    For each root position, slice the sequence from the root position to fw_length for each fw_length allowed. 
    Ex. 
        MIN_LENGTH = 17
        MAX_LENGTH = 22
        pb_start = 607
        1st_root_position = pb_start - 1
        Given the above: 
        primer1 = [590:607] #17 bp primer
        primer2 = [589:607]
        primer3 = [588:607]
        primer4 = [587:607]
        primer5 = [586:607]
        primer6 = [585:607] #22 bp primer

        fw_end --> pb_start - i
        fw_start --> fw_end - fw_len

        {root_pos:pb_start - 1 - i, len:fw_len, seq:[pb_start - i - fw_len :pb_start - i]}
    """
    list_fw_primers = []
    #Search for each root position and each legal primer length
    for i in range(50): 
        for fw_len in range(MIN_PRIMER_LEN, MAX_PRIMER_LEN+1):
            fw_start = pb_start - i - fw_len
            fw_end = pb_start - i
            fw_primer_seq = str(target_seq[fw_start:fw_end])
            if check_primer(fw_primer_seq) is True: 
                list_fw_primers.append(
                    {
                        'root_pos':pb_start - 1 -i,
                        'len':fw_len, 
                        'seq':fw_primer_seq,
                        'tm':float(primer3.calcTm(fw_primer_seq, dv_conc=1.5))
                    }
                )
    return list_fw_primers

def find_rev_primers(target_seq, pb_start, pb_end):
    """
    Function finds all viable REV primers.
    REV primer criteria: 
    1) Between min and max length 
    2) % GC is 30% to 80%
    3) Last five nucleotides at the 3' end contain no more than two G + C residues
    4) No more than 4 consecutive nucleotides within the primer 
    
    Input data: 
    Input data: 
    1) target_seq - str - target sequence - ATCTGATCATGATCATGACTAGTCATGGC
    2) pb_start - int - start index of 5'-end of the probe - 607
    3) pb_end - int - index of 3'-end of the probe (1 past index of last probe nucleotide) - 623 

    Output data: 
    [
        {root_pos:0, len:17, seq:''}
    ]

    Algorithm: 
    The 1st root position is the pb_end index. 
    The last legal root position is the case where: 
        1) fw_primer is at the MIN_PRIMER_LEN
        2) rev_primer is at the MIN_PRIMER_LEN
        3) fw_primer is directly adjacent to probe (no gap)
    To calculate the last legal root position: 
        150 - 2(MIN_PRIMER_LEN) - pb_len
    Let the last legal position be x, min_primer_len be 17, and max_primer len be 22. 
    From 1 to X - (max_primer_len - min_primer_len + 1), there are max_primer_len - min_primer_len + 1 primers. 
    For each subsequent position, the number of primers decreases by 1 until the last position, 
    where the only legal primer length is min_primer_len. 

    {root_pos:pb_end + i, len:rev_len, seq:[(pb_end+i):(pb_end+i+pb_len)]}
    """
    list_rev_primers = []
    len_pb = pb_end - pb_start
    last_root_pos = 150 - 2*MIN_PRIMER_LEN - len_pb
    #Case 1: all primer lengths are available at these root positions
    for i in range(last_root_pos - (MAX_PRIMER_LEN - MIN_PRIMER_LEN + 1)): 
        for rev_len in range(MIN_PRIMER_LEN, MAX_PRIMER_LEN+1): 
            rev_start = pb_end + i
            rev_end = rev_start + rev_len
            rev_primer_seq = str(target_seq[rev_start:rev_end].reverse_complement())
            if check_primer(rev_primer_seq) is True: 
                list_rev_primers.append(
                    {
                        'root_pos':pb_end + i,
                        'len':rev_len,
                        'seq':rev_primer_seq,
                        'tm':float(primer3.calcTm(rev_primer_seq, dv_conc=1.5))
                    }
                )
    #Case 2: end of root positions
    for i in range(MAX_PRIMER_LEN - MIN_PRIMER_LEN + 1):
        rev_start = pb_end + last_root_pos - i
        for rev_len in range(MIN_PRIMER_LEN, MIN_PRIMER_LEN + 1 + i):
            rev_end = rev_start + rev_end
            rev_primer_seq = str(target_seq[rev_start:rev_end].reverse_complement())
            if check_primer(rev_primer_seq) is True: 
                list_rev_primers.append(
                    {
                        'root_pos':pb_end + i,
                        'len':rev_len,
                        'seq':rev_primer_seq,
                        'tm':float(primer3.calcTm(rev_primer_seq, dv_conc=1.5))
                    }
                )
    return list_rev_primers

def blast(seq, blastdb, num_aligned): 
    fasta = tempfile.NamedTemporaryFile(delete=True)
    fasta.write(f">primer\n{str(seq)}".encode())
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

def generate_blast_results_db(list_primers): 
    blast_results = dict()
    for primer in list_primers: 
        primer_id = f"{str(primer['root_pos'])}-{primer['len']}"
        blast_results[primer_id] = blast(primer['seq'], BLASTDB, LEN_DB)
    return blast_results

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

def calculate_primer_pair_sensitivity(fw_blast_results, rev_blast_results, target_accessions): 
    """
    Function calculates the sensitivity of a primer pair given the BLAST results for both. 
    Amplification is defined as an accession where both the forward and reverse primers have 
    100% query coverage and 100% percent identity to the target sequence. 
    Sensitivity is calculated by the following formula: 
        TP / (TP + FN)
        where: 
        TP = succesfully amplified accessions
        FN = target accessions that were not amplified
        TP + FN = total number of target accessions
    """
    #Take only accessions where there was a perfect match --> full query coverage, 100% identity
    fw_perfect_match_results = fw_blast_results.loc[(fw_blast_results['qlen']==fw_blast_results['length'])&(fw_blast_results['pident']==100.0)]
    rev_perfect_match_results = rev_blast_results.loc[(rev_blast_results['qlen']==rev_blast_results['length'])&(rev_blast_results['pident']==100.0)]
    #Retrieve only the accessions list
    fw_perfect_match_accessions = set(fw_perfect_match_results.loc[:,'sacc'])
    rev_perfect_match_accessions = set(rev_perfect_match_results.loc[:,'sacc'])
    amplified_accessions = set(fw_perfect_match_accessions&rev_perfect_match_accessions)
    target_match = 0
    #Count number of target accessions in the amplified accession list 
    for accession in amplified_accessions: 
            if accession in target_accessions: 
                target_match = target_match + 1
    #Calculate sensitivity and return
    sensitivity = target_match/len(target_accessions)
    return sensitivity

def calculate_primer_pair_specificity(fw_blast_results, rev_blast_results, target_accessions, blastdb_len):
    fw_match_accessions = set(fw_blast_results.loc[:,'sacc'])
    rev_match_accessions = set(rev_blast_results.loc[:,'sacc'])
    amplified_accessions = set(fw_match_accessions&rev_match_accessions)
    #Remove every target_accession from the blast_match_accessions list
    for accession in target_accessions:
        if accession in amplified_accessions:  
            amplified_accessions.remove(accession)
    #Calculate specificity
    #Total non-target = all_blast_sequences - target_accessions
    specificity = (blastdb_len - len(amplified_accessions))/(blastdb_len - len(target_accessions))
    return specificity

def find_primer_pairs(list_fw_primers, list_rev_primers, max_tm_diff): 
    """
    Function will pair the fw and rev primers together. 
    The function iterates over the list of FW primers, and pairs them with all viable reverse primers. 
    Viable reverse primers --> a reverse primer such that the amplicon length is equal to or less than 
    150 bp. 

    Algorithm: 
    Amplicon length = (rev_root_pos + rev_len) - fw_root_pos
    Root_pos_limit = FW_root_pos + 150 - MIN_PRIMER_LEN

    """
    def check_tm_diff(fw_tm, rev_tm):
        return abs(fw_tm-rev_tm)

    list_primer_pairs = []

    for fw_primer in list_fw_primers: 
        root_pos_limit = fw_primer['root_pos'] + 150 - MIN_PRIMER_LEN
        for rev_primer in list_rev_primers:
            if (
                rev_primer['root_pos'] < root_pos_limit
                and (rev_primer['root_pos'] + rev_primer['len'] - fw_primer['root_pos']) <= 150
            ):  
                if check_tm_diff(fw_primer['tm'], rev_primer['tm']) <= max_tm_diff:
                    list_primer_pairs.append(
                        {
                            'fw':fw_primer, 
                            'rev':rev_primer,
                            'tm_diff':abs(fw_primer['tm'] - rev_primer['tm']),
                            'average_tm':((fw_primer['tm'] + rev_primer['tm'])/2)
                        }
                )
            else: 
                break
    return list_primer_pairs

def main(target_seq_file_path, target_accessions, pb_start, pb_end, max_tm_diff): 
    target_seq_file = SeqIO.read(target_seq_file_path, 'fasta')
    target_seq = target_seq_file.seq
    list_fw_primers = find_fw_primers(target_seq, pb_start)
    list_rev_primers = find_rev_primers(target_seq, pb_start, pb_end)

    #Generate blast results
    print('BLASTing FW primers...')
    list_fw_blast_results = generate_blast_results_db(list_fw_primers)
    list_rev_blast_results = generate_blast_results_db(list_rev_primers)

    #Generate all legal combination of primer pairs
    list_primer_pairs = find_primer_pairs(list_fw_primers, list_rev_primers, max_tm_diff)

    #Calculate sensitivity and specificity for each primer pair
    for primer_pair in list_primer_pairs: 
        fw_id = f"{str(primer_pair['fw']['root_pos'])}-{primer_pair['fw']['len']}"
        rev_id = f"{str(primer_pair['rev']['root_pos'])}-{primer_pair['rev']['len']}"
        fw_blast_result = list_fw_blast_results[fw_id]
        rev_blast_result = list_rev_blast_results[rev_id]
        primer_pair['sensitivity'] = calculate_primer_pair_sensitivity(fw_blast_result, rev_blast_result, target_accessions)
        primer_pair['specificity'] = calculate_primer_pair_specificity(fw_blast_result, rev_blast_result, target_accessions, LEN_DB)
        primer_pair['score'] = primer_pair['sensitivity'] + primer_pair['specificity']

    #Sort
    list_primer_pairs.sort(key=lambda x: (x['score'], x['average_tm'], x['tm_diff']), reverse=True)

    #Outputs: 
    #1) blast_csv for each blast result for each primer
    #2) file listing all of the primer pairs

    blast_folder_path = target_seq_file_path.with_name('primer_blast')
    blast_folder_path.mkdir()
    for blast_result_key in list_fw_blast_results: 
        blast_output_path = blast_folder_path.joinpath(f"{blast_result_key}_fw.csv")
        blast_output_file = open(blast_output_path, 'w')
        list_fw_blast_results[blast_result_key].to_csv(blast_output_file)
        blast_output_file.close()
    for blast_result_key in list_rev_blast_results: 
        blast_output_path = blast_folder_path.joinpath(f"{blast_result_key}_rev.csv")
        blast_output_file = open(blast_output_path, 'w')
        list_rev_blast_results[blast_result_key].to_csv(blast_output_file)
        blast_output_file.close()

    #CSV primer pair output
    primer_pair_data = []
    for primer_pair in list_primer_pairs: 
        primer_pair_data.append(
            (
                primer_pair['fw']['root_pos'],
                primer_pair['fw']['seq'],
                primer_pair['fw']['len'],
                primer_pair['fw']['tm'],
                primer_pair['rev']['root_pos'],
                primer_pair['rev']['seq'],
                primer_pair['rev']['len'],
                primer_pair['rev']['tm'],
                primer_pair['sensitivity'],
                primer_pair['specificity'],
                primer_pair['score'],
            )
        )
    csv_output_path = target_seq_file_path.with_name('primer_pairs.csv')
    csvfile = open(csv_output_path, 'w', newline='')
    csvwriter = csv.writer(csvfile)
    csvwriter.writerow(
        (
            'fw_root',
            'fw_len',
            'fw_seq',
            'fw_tm',
            'rev_root',
            'rev_len',
            'rev_seq',
            'rev_tm',
            'sens',
            'spec',
            'score'
        )
    )
    csvwriter.writerows(primer_pair_data)
    csvfile.close()

if __name__ == '__main__':
    #target_seq_file_path = 'ssu_ppenetrans.fasta'
    #pb_start = 670
    #pb_end = 689
    #input_file = open('list_accessions.txt', 'r')
    target_seq_file_path, target_accessions_path, pb_start, pb_end, tm_diff = parse_args()
    target_accessions = []
    input_file = open(target_accessions_path, 'r')
    for line in input_file: 
        target_accessions.append(line.strip('\n'))
    input_file.close()
    main(target_seq_file_path, target_accessions, pb_start, pb_end, tm_diff)