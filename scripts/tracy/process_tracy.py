import json
import argparse
import pathlib
import csv

def parse_metadata(all_seq_data): 
    trim_parameters = []
    for cow in all_seq_data:
        seq_len = len(all_seq_data[cow]['primarySeq'])
        trimleft = all_seq_data[cow]['meta']['arguments']['trimLeft']
        trimright = all_seq_data[cow]['meta']['arguments']['trimRight']
        adjusted_seq_len = seq_len - trimleft - trimright
        trim_parameters.append((cow, seq_len, adjusted_seq_len, trimleft, trimright))
    return trim_parameters

def parse_variants(all_seq_data): 
    variants = []
    for cow in all_seq_data: 
        list_variants = all_seq_data[cow]['variants']['rows']
        if len(list_variants) > 0: 
            for variant in list_variants:
                variant.insert(0, cow)
                variants.append(variant)
    return variants

def output_csv(data, file_name, file_path, header): 
    output_path = file_path.with_name(file_name+'.csv')
    output_file = open(output_path, 'w', newline='')
    csv_writer = csv.writer(output_file, dialect='excel')
    csv_writer.writerow(header)
    csv_writer.writerows(data)
    output_file.close()

def parse_args(): 
    parser = argparse.ArgumentParser(description='Extract variant calls from Tracy JSON file')
    parser.add_argument(
        'json_path',
        type=pathlib.Path,
        action='store',
        help='Path to .json files'
    )
    args = parser.parse_args()
    return args.json_path

def main(json_path): 
    #Get all json files in the dir
    all_seq_data = dict()
    if json_path.is_dir() is True: 
        for json_file_path in json_path.glob('*.json'): 
            cow = json_file_path.stem.split('_')[1]
            input_file = open(json_file_path, 'r')
            seq_data = json.load(input_file)
            all_seq_data[cow] = seq_data
            input_file.close()
    
    metadata = parse_metadata(all_seq_data)
    variants = parse_variants(all_seq_data)

    #Output
    metadata_header = [
        'Sample',
        'Raw Sequence Length',
        'Trimmed Sequence Length', 
        'Left Trim', 
        'Right Trim',
    ]
    variants_header = [
        'Sample',
        'chr',
        'pos',
        'id',
        'ref',
        'alt',
        'qual',
        'filter',
        'type',
        'genotype',
        'basepos',
        'signalpos',
    ]
    output_csv(metadata, 'metadata', json_path, metadata_header)
    output_csv(variants, 'variants', json_path, variants_header)

if __name__ == '__main__': 
    json_path = parse_args()
    main(json_path)