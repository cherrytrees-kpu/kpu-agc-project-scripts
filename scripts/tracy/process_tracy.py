import json
import argparse
import pathlib

def parse_metadata(all_seq_data): 
    trim_parameters = []
    for cow in all_seq_data:
        seq_len = len(all_seq_data[cow]['primarySeq'])
        trimleft = all_seq_data[cow]['meta']['arguments']['trimLeft']
        trimright = all_seq_data[cow]['meta']['arguments']['trimRight']
        trim_parameters.append((cow, trimleft, trimright))
    print(trim_parameters)

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
    
    parse_metadata(all_seq_data)

if __name__ == '__main__': 
    json_path = parse_args()
    main(json_path)