import csv
import argparse
import pathlib

def parse_args():
    parser = argparse.ArgumentParser(description='Prepare a file for the Phusion Tm Calculator')
    parser.add_argument(
        'primer_path',
        action='store',
        type=pathlib.Path, 
        help='Path to .tab file containing the primers, fw primer sequence, rev primer sequence' 
    )
    args = parser.parse_args()
    return (args.primer_path)

def main(primer_path):
    degen_nucleotide = {
        'Y':['C', 'T'],
        'R':['A', 'G']
    }
    primer_data = []
    csv_file = open(primer_path, 'r')
    csv_reader = csv.reader(csv_file, dialect='excel', delimiter='\t')
    for line in csv_reader: 
        primer_id = line[0]
        fw_primer = line[1]
        rev_primer = line[2]
        primer_data.append((primer_id, fw_primer, rev_primer))
    output_path = primer_path.with_name('phusion_tm_file.txt')
    output_file = open(output_path, 'w')
    for primer in primer_data: 
        output_file.write(
            f'{primer[0]}-FW {primer[1]} ; {primer[0]}-REV {primer[2]}\n'
        )
    output_file.close()

if __name__ == '__main__':
    primer_path = parse_args()
    main(primer_path)