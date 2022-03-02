#!/usr/bin/env python3
"""
Author : Erick Samera
Date   : 2022-02-05
Purpose: Parse tracy JSON files and produce summary .xlsx sheet.
"""

import json
import pathlib
import argparse
import random
import time
from typing import NamedTuple
import pandas as pd

class Args(NamedTuple):
    """ Command-line arguments """
    input_path: pathlib.Path
    make_template: bool
    use_template: pathlib.Path

    verbose: int
    output_csv: bool

# --------------------------------------------------
def get_args() -> Args:
    """ Get command-line arguments """

    parser = argparse.ArgumentParser(
        description='Parse tracy JSON files and produce a summary .xlsx spreadsheet.',
        usage='%(prog)s [--make_template/--summary] [options] [path]',
        formatter_class=argparse.RawTextHelpFormatter
        )

    parser.add_argument(
        'input_path',
        metavar = 'input_path',
        type = pathlib.Path,
        help = '\n'.join([
            'path to perform function:',
            '--make_template    output for template',
            '--template_input   dir containing dirs of tracy JSON files'
            ])
        )

    main_arg = parser.add_mutually_exclusive_group(required=True)

    main_arg.add_argument(
        '--make_template',
        action='store_true',
        help='create a template and exit'
        )

    main_arg.add_argument(
        '--use_template',
        metavar='FILE',
        type=pathlib.Path,
        help='summarizes data using specified template'
        )

    group_debug = parser.add_argument_group('debugging')

    group_debug.add_argument(
        '-v',
        '--verbose',
        action='count',
        default=0,
        help='set verbosity',
        )

    group_debug.add_argument(
        '--output_csv',
        action='store_true',
        help='output .csv for each JSON group',
        )

    args = parser.parse_args()

    return Args(args.input_path, args.make_template, args.use_template, args.verbose, args.output_csv)

# --------------------------------------------------
def main() -> None:
    """ Do the thing """

    # get command-line args
    args = get_args()
    input_path_arg = args.input_path.resolve()
    make_template_arg = args.make_template
    use_template_arg = args.use_template.resolve() if args.use_template else None

    verbose_arg = args.verbose
    output_csv_arg = args.output_csv

    # define output_path
    output_path = input_path_arg

    # create a template if the template flag was specified
    if make_template_arg:
        print_runtime('created template file')
        print('refer to the template below for formatting: \n')
        generate_template_file(output_path)
        print('\naliases not case-sensitive \n')

    # carry on with the rest of tracy parse
    elif use_template_arg:
        region_dir_genotypes = {}

        # for each region, generate a set of genotypes
        for region_dir in input_path_arg.iterdir():
            if region_dir.is_dir():
                region_dir_genotypes[region_dir.stem] = generate_genotype_DataFrame(use_template_arg, region_dir)
                print_runtime(f'Generated genotypes DataFrame for {region_dir.name}')
                if output_csv_arg:
                    generate_genotype_DataFrame(use_template_arg, region_dir).to_csv(output_path.joinpath(f'{region_dir.name}_genotypes.csv'), index=False)
                    print_runtime(f'Created file {region_dir.name}_genotypes.csv')

        write_genotype_DataFrames_to_excel(region_dir_genotypes, output_path)
        print_runtime(f'Finished writing excel file!')

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
        f_neutral_underline = workbook.add_format({
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
        f_bad_underline = workbook.add_format({
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
        f_sample_name_column = workbook.add_format({
            'align' : 'center',
            'valign' : 'vcenter',
            'bold' : True,
            'bottom' : 1 })
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

            list_of_genotype_DataFrames[gene].to_excel(writer, sheet_name=gene, index=False)
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

            # change the column formatting sizes and hide some
            worksheet.set_column(0, 1, 14, f_align_left)
            worksheet.set_column(1, 2, 12, f_align_left)
            worksheet.set_column(3, len(list_of_genotype_DataFrames[gene].columns.tolist()) + 1, 14.3, f_align_center)
            worksheet.set_column(1, 1, options={'hidden' : True})

            # add header formatting
            for col_num, col_value in enumerate(list_of_genotype_DataFrames[gene].columns.values):
                worksheet.write(0, col_num, col_value, f_header)					                                # header
                worksheet.write(1, col_num, list_of_genotype_DataFrames[gene].iat[0, col_num], f_ref)		        # subheader

            # add sample_name formatting
            for sample_i, sample_name in enumerate(list_of_genotype_DataFrames[gene]['sample_name']):
                worksheet.write(sample_i+1, 0, sample_name, f_sample_name_column)

            # iterate through the rows and do actions
            for row_i, row_val in genetype_DataFrame_for_excel[::3].iterrows():
                # hide the second and third row after the reference
                worksheet.set_row(row_i+2, None, options={'hidden' : True})
                worksheet.set_row(row_i+3, None, options={'hidden' : True})

                # for each column, get the column reference
                for col_i, col_v in enumerate(list_of_genotype_DataFrames[gene].columns.values):
                    if col_i > 2:

                        reference_genotype = str(genetype_DataFrame_for_excel.at[0, col_v]).strip()
                        individual_genotype = str(row_val[col_v]).strip()
                        #print(reference_genotype, individual_genotype)

                        middle_line = 12
                        bottom_line = 24

                        if row_i == middle_line:
                            worksheet.write(row_i + 1, col_i, genetype_DataFrame_for_excel.at[row_i, col_v], f_bos_divider)
                        if row_i == bottom_line:
                            worksheet.write(row_i + 1, col_i, genetype_DataFrame_for_excel.at[row_i, col_v], f_bos_divider)

                        f_neutral_type = f_neutral_underline if (row_i == middle_line or row_i == bottom_line) else f_neutral
                        f_bad_type = f_bad_underline if (row_i == middle_line or row_i == bottom_line) else f_bad

                        # if the genotype isn't the reference and isn't a failed sequence, highlight this cell as yellow
                        if individual_genotype not in (reference_genotype, '*'):
                            worksheet.write(row_i + 1, col_i, genetype_DataFrame_for_excel.at[row_i, col_v], f_neutral_type)

                        # if the genotype is a failed sequence, highlight this cell as red
                        elif individual_genotype == '*':
                            worksheet.write(row_i + 1, col_i, '*', f_bad_type)

            worksheet.merge_range('C5:C14', 'Bos taurus', f_species_merge)
            worksheet.merge_range('C17:C26', 'Bos indicus', f_species_merge)

            print_runtime(f'Created worksheet for {gene}')

def generate_genotype_DataFrame(input_template, input_region_dir) -> pd.DataFrame:
    """
    Function will parse_JSON files in a given directory and create a dict containing genotypes of each sample against a reference
    """

    template_dict = read_template(input_template)
    last_orig_column_pos = len(template_dict[0].keys()) - 1

    for json_file in input_region_dir.glob('*.json'):

        # update the reference genotypes
        for SNP in parse_json(json_file)['reference_SNP']:
            template_dict[0].update(SNP)
            template_dict[0]['seq'] = True

        # update the sample genotypes
        sample_i = validate_sample_name(json_file.stem.split('_')[1], template_dict)['name_position']
        for SNP in parse_json(json_file)['sample_SNP']:
            template_dict[sample_i].update(SNP)

        # update seq flag since tracy json exists
        template_dict[sample_i]['seq'] = True

    # clean up the DataFrame before returning it
    template_DataFrame = pd.DataFrame(pd.DataFrame.from_dict(template_dict, orient='index'))
    template_DataFrame.drop(columns='aliases', inplace=True)
    sorted_columns = list(template_DataFrame.columns[:last_orig_column_pos]) + sorted(template_DataFrame.columns[last_orig_column_pos:])
    template_DataFrame = template_DataFrame[sorted_columns]

    # replace all blanks before exploding
    # iterate through the tracy variant-called columns
    # define the reference genotype and compare the sample genotype to that
    # also check whether the sequence failed
    for column in list(template_DataFrame.columns)[last_orig_column_pos:]:
        reference_genotype = template_DataFrame.at[0, column][0]
        for row_num, _ in template_DataFrame.iterrows():
            sample_genotype = template_DataFrame.at[row_num, column]
            if not isinstance(sample_genotype, list):
                if template_DataFrame.at[row_num, 'seq'] and pd.isna(sample_genotype):
                    template_DataFrame.at[row_num, column] = [reference_genotype, None, None]
                elif not template_DataFrame.at[row_num, 'seq'] and pd.isna(sample_genotype):
                    template_DataFrame.at[row_num, column] = ['*', None, None]

    try:
        # explode
        template_DataFrame = template_DataFrame.explode(list(template_DataFrame.columns)[last_orig_column_pos:])
    except ValueError:
        print(template_DataFrame)

    return template_DataFrame

def parse_json(query_json) -> dict:
    """
    Function will parse json files produced by tracy and return a list of sample genotypes and corresponding reference genotypes
    """

    with open(f'{query_json}', 'r', encoding='UTF-8') as tracy_json:

        # initialize the json file from tracy
        data = json.load(tracy_json)

        # make a dataframe of the variants
        variant_data = pd.DataFrame(data['variants']['rows'], columns=data['variants']['columns'])

        # initialize the genotype and reference for each variant position
        sample_genotypes_list = []
        reference_genotypes_list = []

        for _, row_val in variant_data.iterrows() if not variant_data.empty else []:
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
                        row_val['ref'],
                        None,
                        None,
                        ]})

    return {'sample_SNP': sample_genotypes_list, 'reference_SNP': reference_genotypes_list}

def validate_sample_name(input_name, template_dict) -> dict:
    """
    Function will validate the input name with a sample list of defined names and aliases
    """

    verified_sample_name = [key for key, value in template_dict.items() if input_name.lower() in value['aliases']]

    if verified_sample_name:
        return {'name_value': template_dict[verified_sample_name[0]]['sample_name'], 'name_position': verified_sample_name[0]}

def generate_template_file(output_path) -> None:
    """
    Function will generate a .tab template file for use with the --template_input flag.
    """

    time_string = time.strftime("%Y_%m_%d_%H%M%S", time.localtime(time.time()))

    with open(output_path.joinpath(f'template_{time_string}.tab'), 'w', encoding='UTF-8') as template_file:
        headers = [
            'sample_name',
            'aliases'
        ]
        template_file.write('\t'.join(headers))

    example_data = {
        0 : [
            ('Erick', ['erick']),
            ('RJ', ['robert john', 'rj']),
            ('Michael', ['cherry_trees']),
            ('Monique', ['goof']),
            ('Daleena', ['williams'])],
        1 : [
            ('AliBarber', ['alibarber']),
            ('Cochise', ['cochise']),
            ('Sansao', ['sansao']),
            ('Slugger', ['slugger', 'slogger']),
            ('LLNKune', ['llnkune', 'llnkure']),
            ('Nagaki', ['nagaki']),
            ('Raider', ['raider']),
            ('Rendition', ['rendition', 'renditium']),
        ]}

    random_example_data = random.choices(example_data)[0]
    example_headers = ['sample_name', 'aliases']

    example_DataFrame = pd.DataFrame(random_example_data, columns=example_headers)
    example_DataFrame['param+1'] = ''
    print(example_DataFrame)

def read_template(input_template_file) -> dict:
    """
    Function will read template file and output a dict that will be easier to deal with.
    """

    template_file_DataFrame = pd.DataFrame(pd.read_csv(input_template_file, delimiter='\t'))

    # create a reference row at the top of the DataFrame
    template_file_DataFrame.loc[-1] = ['#reference',] + (['#ref',] * (len(template_file_DataFrame.columns) - 1))
    template_file_DataFrame.index = template_file_DataFrame.index + 1
    template_file_DataFrame = template_file_DataFrame.sort_index()

    # fix the aliases column
    for alias_i, alias_value in enumerate(template_file_DataFrame['aliases']):
        template_file_DataFrame.at[alias_i, 'aliases'] = [alias.lower().strip() for alias in alias_value[1:-1].split(',')]

    # insert sequence flag for successful sequencing
    template_file_DataFrame.insert(1, 'seq', [False,] * len(template_file_DataFrame))

    return template_file_DataFrame.to_dict(orient='index')

def print_runtime(action) -> None:
    """ Print the time and some defined action. """
    print(f'[{time.strftime("%Y-%m-%d %H:%M:%S", time.localtime())}] {action}')

# --------------------------------------------------
if __name__ == '__main__':
    main()