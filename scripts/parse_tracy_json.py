#!/usr/bin/env python3
"""
Author : Erick Samera
Date   : 2022-02-03
Purpose: Parse tracy JSON files and produce summary .xlsx sheet.
"""

import argparse
from typing import NamedTuple
import json, pathlib, time
import pandas as pd

class Args(NamedTuple):
    """ Command-line arguments """
    json_file_path: pathlib.Path
    output_dir: bool
    input_type: bool
    csv: bool

# --------------------------------------------------
def get_args() -> Args:
    """ Get command-line arguments """

    parser = argparse.ArgumentParser(
        description='Parse tracy JSONs and produce a summary excel sheet',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument('json_file_path',
                        metavar='json_file_path',
                        type=pathlib.Path,
                        help='Path of directory containing tracy JSONs (or directories of gene JSONs)')

    parser.add_argument('-o',
                        '--output_dir',
                        help="flag whether directory 'output' will be created",
                        action='store_false')

    parser.add_argument('-i',
                        '--input_type',
                        help="flag whether Michael put genes into individual folders",
                        action='store_false')
    
    parser.add_argument('--csv',
                        help="flag whether .csv files will be produced for each gene (mostly for debugging)",
                        action='store_true')

    args = parser.parse_args()

    return Args(args.json_file_path, args.output_dir, args.input_type, args.csv)

# --------------------------------------------------
def main() -> None:

    start_main = time.time()

    # define args
    args = get_args()
    json_file_path_arg = args.json_file_path.resolve()
    output_dir_arg = args.output_dir
    input_type_arg = args.input_type
    csv_arg = args.csv

    # check output_dir flag, if true: make an output dir
    if output_dir_arg:
        output_flag = 'tracy_json_parse'
        output_path = json_file_path_arg.joinpath(f"{output_flag}_output-{time.strftime('%Y_%m_%d-%H%M%S', time.localtime(time.time()))}")
        output_path.mkdir(parents=True, exist_ok=True)
    else:
        output_path = json_file_path_arg

    sample_list = {
        '#reference': {
                'aliases': ['reference', 'ref'],
                'species': '#ref'
                },
        'AliBarber': {
                'aliases': ['alibarber'],
                'species': 'Bos taurus'
                },
        'Cochise': {
                'aliases': ['cochise'],
                'species': 'Bos taurus'
                },
        'Sansao': {
                'aliases': ['sansao'],
                'species': 'Bos taurus'
                },
        'Slugger': {
                'aliases': ['slugger', 'slogger'],
                'species': 'Bos taurus'
                },
        'LLNKune': {
                'aliases': ['llnkune', 'llnkure'],
                'species': 'Bos indicus'
                },
        'Nagaki': {
                'aliases': ['nagaki'],
                'species': 'Bos indicus'
                },
        'Raider': {
                'aliases': ['raider'], 
                'species': 'Bos indicus'
                },
        'Renditium': {
                'aliases': ['renditium', 'rendition'],
                'species': 'Bos indicus'
                }
        }

    list_of_genotype_DataFrames = {}

    # check input_type flag, if true: 
    if input_type_arg:
        for gene_dir in json_file_path_arg.iterdir():
            if gene_dir.is_dir():
                for sample_json in gene_dir.glob('*.json'):
                    query_path = gene_dir
                    list_of_genotype_DataFrames.update({gene_dir.stem: generate_genotype_DataFrame(sample_list, gene_dir.stem, query_path, output_path, csv_arg)})
    else:
        list_of_genes = list(set([gene.stem.split('_')[0] for gene in json_file_path_arg.glob('*.json')]))
        for gene in list_of_genes:
            query_path = json_file_path_arg
            list_of_genotype_DataFrames.update({gene: generate_genotype_DataFrame(sample_list, gene, query_path, output_path, csv_arg)})
    write_genotype_DataFrames_to_excel(list_of_genotype_DataFrames, output_path)
    print_runtime(f'Produced summary tracy genotyping excel file.')

def generate_genotype_DataFrame(sample_list, gene_ID, query_path, output_path, csv_arg):
    """
    Function will generate a genotype DataFrame with genotyping template
    """
    SNP_data = generate_template(sample_list)
            
    for sample_json in query_path.glob(f'{gene_ID}*.json'):
        
        gene = sample_json.stem.split('_')[0]
        sample = validate_sample_name(sample_json.stem.split('_')[1], sample_list)

        results = parse_json(sample_json, gene, sample)
        SNP_data[sample]['seq'] = True
        SNP_data['#reference']['seq'] = True

        for i in results[0]:
            SNP_data[sample].update(i)
        
        for i in results[1]:
            SNP_data['#reference'].update(i)

    SNP_DataFrame = pd.DataFrame.from_dict(SNP_data, orient='index')

    for column in SNP_DataFrame.columns[2:]:
        reference_genotype = SNP_DataFrame.at['#reference', column][0]
        
        for row_num, row_val in SNP_DataFrame[1:].iterrows():
            if not isinstance(SNP_DataFrame.at[row_num, column], list) and SNP_DataFrame.at[row_num, 'seq']:
                SNP_DataFrame.at[row_num, column] = [reference_genotype, None, None]
            elif not SNP_DataFrame.at[row_num, 'seq']:
                SNP_DataFrame.at[row_num, column] = ['*', None, None]

    if csv_arg:
        SNP_DataFrame.to_csv(output_path.joinpath(f'{gene}.csv'))
        print_runtime(f'Produced {gene}.csv')
    
    if SNP_DataFrame.columns.tolist()[2:]:
        SNP_DataFrame = SNP_DataFrame.explode(SNP_DataFrame.columns.tolist()[2:])
    return SNP_DataFrame

def write_genotype_DataFrames_to_excel(list_of_genotype_DataFrames, output_path) -> None:
    """
    Function will write each of the genotypes into a summary excel sheet with each gene in their own sheet
    """

    with pd.ExcelWriter(output_path.joinpath('#summary_tracy_results.xlsx'), engine='xlsxwriter') as writer:
        workbook = writer.book

        # excel formatting styles
        # --------------------------------------------------
        f_align_left = workbook.add_format({
            'align': 'left',
            'valign': 'vcenter',})
        f_align_center = workbook.add_format({
            'align': 'center',
            'valign': 'vcenter',})
        f_neutral = workbook.add_format({
            'align' : 'center',
            'valign' : 'vcenter',
            'bg_color' : '#ffeb9c',
            'font_color' : '#9c5700'})
        f_neutral_under = workbook.add_format({
            'align' : 'center',
            'valign' : 'vcenter',
            'bg_color' : '#ffeb9c',
            'font_color' : '#9c5700',
            'bottom' : 1 })
        f_bad = workbook.add_format({
            'align' : 'center',
            'valign' : 'vcenter',
            'bg_color' : '#ffc7ce',
            'font_color' : '#9c0006'})
        f_bad_under = workbook.add_format({
            'align' : 'center',
            'valign' : 'vcenter',
            'bg_color' : '#ffc7ce',
            'font_color' : '#9c0006',
            'bottom' : 1 })
        f_header = workbook.add_format({
            'align' : 'center',
            'valign' : 'vcenter',
            'bold' : True,
            'bottom' : 6 })
        f_ref = workbook.add_format({
            'align' : 'center',
            'valign' : 'vcenter',
            'bold' : True,
            'bottom' : 1 })
        f_bos_divider = workbook.add_format({
            'align' : 'center',
            'valign' : 'vcenter',
            'bottom' : 1 })
        f_species_merge = workbook.add_format({
            'bold': 1,
            'border': 1,
            'align': 'center',
            'valign': 'vcenter',
            'rotation': 90})

        # for each each gene in the list of genotype DataFrames, set up an excel sheet
        for gene in list_of_genotype_DataFrames:

            # deal with blank tables because tracy didn't find any variants
            if len(list_of_genotype_DataFrames[gene]) == 9:
                list_of_genotype_DataFrames[gene][''] = [None]*9
                list_of_genotype_DataFrames[gene].at['#reference', ''] = [None]
                for row_i, row_val in list_of_genotype_DataFrames[gene].iterrows():
                    if row_i != '#reference':
                        list_of_genotype_DataFrames[gene].at[row_i,''] = [None, None, None]
                list_of_genotype_DataFrames[gene] = list_of_genotype_DataFrames[gene].explode('')

            list_of_genotype_DataFrames[gene].to_excel(writer, sheet_name=gene, index = True)
            genetype_DataFrame_for_excel  = list_of_genotype_DataFrames[gene].reset_index()

            paternal_expressed = ['GNAS', 'MEST', 'NNAT', 'PEG10', 'PEG3', 'PLAGL1', 'RTL1', 'SNRPN', 'XIST']
            maternal_expressed = ['H19', 'IGF2R']
            # make formatting changes 
            # --------------------------------------------------
            worksheet = writer.sheets[gene]

            # colour sheet depending on whether it's maternal or paternal expressed
            if gene.split('-')[0] in paternal_expressed:
                worksheet.set_tab_color('#89CFF0')
            elif gene.split('-')[0] in maternal_expressed:
                worksheet.set_tab_color('#F2ACB9')

            # change the column formatting
            worksheet.set_column(0, 1, 10, f_align_left)
            worksheet.set_column(1, 2, 7, f_align_left)
            worksheet.set_column(3, len(list_of_genotype_DataFrames[gene].columns.tolist()) + 1, 14.3, f_align_center)
            worksheet.set_column(2, 3, options={'hidden' : True})

            # write the header
            for col_num, col_value in enumerate(list_of_genotype_DataFrames[gene].columns.values):
                worksheet.write(0, col_num + 1, col_value, f_header)					                                # header
                worksheet.write(1, col_num + 1, list_of_genotype_DataFrames[gene].iat[0, col_num], f_ref)		        # subheader
            
            # iterate through the rows and do actions
            for row_i, row_val in genetype_DataFrame_for_excel[1::3].iterrows():
                # hide the second and third row after the reference
                worksheet.set_row(row_i+2, None, options={'hidden' : True}) 
                worksheet.set_row(row_i+3, None, options={'hidden' : True})

                # for each column, get the column reference
                for col_i, col_v in enumerate(list_of_genotype_DataFrames[gene].columns.values):
                    if col_i > 1:

                        reference_genotype = str(genetype_DataFrame_for_excel.at[0, col_v]).strip()
                        individual_genotype = str(row_val[col_v]).strip()

                        if row_i == 10: 
                            worksheet.write(row_i + 1, col_i + 1, genetype_DataFrame_for_excel.at[row_i, col_v], f_bos_divider)
                        if row_i == 22: 
                            worksheet.write(row_i + 1, col_i + 1, genetype_DataFrame_for_excel.at[row_i, col_v], f_bos_divider)

                        f_neutral_type = f_neutral_under if (row_i == 10 or row_i == 22) else f_neutral
                        f_bad_type = f_bad_under if (row_i == 10 or row_i == 22) else f_bad

                        # if the genotype isn't the reference and isn't a failed sequence, highlight this cell as yellow
                        if individual_genotype != reference_genotype and individual_genotype != '*':
                            worksheet.write(row_i + 1, col_i + 1, genetype_DataFrame_for_excel.at[row_i, col_v], f_neutral_type)

                        # if the genotype is a failed sequence, highlight this cell as red
                        elif individual_genotype == '*':
                            worksheet.write(row_i + 1, col_i + 1, '*', f_bad_type)

            worksheet.merge_range('B3:B12', 'Bos taurus', f_species_merge)
            worksheet.merge_range('B15:B24', 'Bos indicus', f_species_merge)

def parse_json(query_ab1, gene, name):
    """
    Function will parse json files produced by tracy and return a list of sample genotypes and corresponding reference genotypes
    """

    with open(f'{query_ab1}') as tracy_json:

        # initialize the json file from tracy
        data = json.load(tracy_json)

        # make a dataframe of the variants
        variant_data = pd.DataFrame(data['variants']['rows'], columns=data['variants']['columns'])

        # initialize the genotype and reference for each variant position
        sample_genotypes_list = []
        reference_genotypes_list = []
        
        for row_i, row_val in variant_data.iterrows() if not variant_data.empty else []:
            if row_val['filter'] == 'PASS':
                
                # define the reference position relative to the reference genome
                genomic_pos = f"{row_val['chr']}:{row_val['pos']}"

                # pretty up the genotype information
                if row_val['genotype'] == 'het.':
                    genotype = f"{row_val['ref']}/{row_val['alt']}"

                elif row_val['genotype'] == 'hom. ALT':
                    genotype = f"{row_val['alt']}"

                # update the genotype for the position
                sample_genotypes_list.append({
                    genomic_pos: [
                        genotype,
                        int(row_val['basepos']),
                        row_val['genotype']
                        ]})
                
                reference_genotypes_list.append({
                    genomic_pos: [
                        row_val['ref']
                        ]})

    return sample_genotypes_list, reference_genotypes_list

def validate_sample_name(query_name: str, sample_list):
    """
    Function will validate the query name with a sample list of defined names and aliases
    """

    verified_sample_name = [key for key, value in sample_list.items() if query_name.lower() in value['aliases']]

    if verified_sample_name:
        return verified_sample_name[0]
    else:
        return None

def generate_template(sample_list):
    """
    Function will generate a template of the genotype table
    """

    genotypes_template = {}

    for sample in sample_list:
        genotypes_template.update({
            sample : {'species': sample_list[sample]['species'], 'seq': False}
            })

    return genotypes_template

def print_runtime(action):
    print(f'[{time.strftime("%Y-%m-%d %H:%M:%S", time.localtime())}] {action}')
# --------------------------------------------------
if __name__ == '__main__':
    main()
