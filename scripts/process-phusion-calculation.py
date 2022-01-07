import csv
import argparse
import pathlib

def parse_args():
    parser = argparse.ArgumentParser(description='Process csv output from Thermofisher Tm calculator')
    parser.add_argument(
        'csv_path',
        action='store',
        type=pathlib.Path, 
        help='Path to .csv file containing primer information' 
    )
    args = parser.parse_args()
    return (args.csv_path)

def main(csv_path):
    csv_file = open(csv_path, 'r')
    csv_reader = csv.reader(csv_file)
    data = []

    for line in csv_reader: 
        if line[0] != '': 
            data.append(line)
    
    output_path = csv_path.with_name('tm_calcs.csv')
    output_file = open(output_path, 'w', newline='')

    csv_writer = csv.writer(output_file)
    csv_writer.writerows(data)
    output_file.close()

if __name__ == '__main__':
    csv_path = parse_args()
    main(csv_path)