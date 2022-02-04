#!/usr/bin/env python3
"""
Author : Erick Samera
Date   : 2021-12-30 -> 2022-02-03
Purpose: Performs BiSearch primer design for all .fasta files in a given folder.
"""

import argparse
from typing import NamedTuple

import pathlib
import re
import random
import time
import pandas as pd
from Bio import SeqIO

try:
    import mechanicalsoup
except ImportError as e:
    print('This relies on MechanicalSoup')

class Args(NamedTuple):
    """ Command-line arguments """

    fasta_file_path: pathlib.Path
    output_dir: bool

    bis: bool
    strand: str
    max: int

    sub: int
    tries: int

    verbose: bool
    seed: int

# --------------------------------------------------
def get_args() -> Args:
    """ Get command-line arguments """

    parser = argparse.ArgumentParser(
        description='Program parses .fasta files and generates primers for bisulfite-converted template via BiSearch.',
        usage='%(prog)s [options] [fasta_file_path]',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument('fasta_file_path',
        metavar='fasta_file_path',
        type=pathlib.Path,
        help='.fasta file or directory containing .fasta files')

    parser.add_argument('-o',
        '--output_dir',
        help="flag whether directory 'output' will be created",
        action='store_true')

    group_BiSearch = parser.add_argument_group('BiSearch parameters')

    group_BiSearch.add_argument('--bis',
        help="bisulfite-converted",
        action='store_true')

    group_BiSearch.add_argument('--strand',
        help="specify target strand after bisulfite-conversion {sense= 'p', antisense = 'm'}",
        choices=['p', 'm'],
        type=str,
        default=None)

    group_BiSearch.add_argument('--max',
        help="max length for PCR",
        metavar='<n>',
        type=int,
        default=599)

    group_subsequence = parser.add_argument_group('subsequence search options')

    group_subsequence.add_argument('-s',
        '--sub',
        help='divide the sequence into <n> subsequences to pull tries from',
        metavar='<n>',
        type=int,
        default=1)

    group_subsequence.add_argument('-t',
        '--tries',
        help='number of tries to attempt in each subsequences',
        metavar='<n>',
        type=int,
        default=1)

    group_subsequence.add_argument('-f',
        '--focus',
        help='focus search only on CpG-dense region',
        action='store_false')

    group_debug = parser.add_argument_group('debug')

    group_debug.add_argument('-v',
        dest='verbose',
        action='count',
        default=0,
        help='verbose output for debugging')

    group_debug.add_argument('--seed',
        help='seed random',
        metavar='<n>',
        type=int,
        default=990318)

    args = parser.parse_args()

    if args.strand and not args.bis:
        parser.error('--strand is only applicable if bisulfite-converted.')

    return Args(args.fasta_file_path, args.output_dir, args.bis, args.strand, args.max, args.sub, args.tries, args.verbose, args.seed)

# --------------------------------------------------
class BiSearch_condition:
    """
    BiSearch conditions to be queried with a browser instance
    """
    BiSearch_form_mapping = (
            'seq',		# query sequence
            'bis',		# bisulfite status (True = bisulfite-converted)
            'sens',		# strand (only available if bisulfite-converted)
            'fbeg',		# forward primer, beginning of primer query
            'fend',		# forward primer, end of primer query
            'rbeg',		# reverse primer, beginning of primer query
            'rend',		# reverse primer, end of primer query
            'len',		# maximum PCR length

            'pc',		# primer concentration (umol)
            'kc',		# potassium concentration (mmol)
            'mgc',		# magnesium concentration (mmol)

            'nstore',	# results list size (still sorted by score, not necessary)

            'db')		# database for ePCR

    def __init__(
            self,
            sequence: str = '',             # query sequence to search for primers
            bisulfite: bool = False,        # bisulfite conversion status, True = bisulfited converted
            strand = None,                  # only applicable after bislfite conversion, 'p' = (+) strand, 'm' = (-) strand
            forward_start: int = None,      # forward primer, start of query region
            forward_end: int = None,        # forward primer, end of query region
            reverse_start: int = None,      # reverse primer, start of query region
            reverse_end: int = None,        # reverse primer, end of query region
            PCR_max_length: int = 599,      # maximum length (nucleotides) of PCR amplicon to query

            primer_conc: float = 2,         # concentration of primers (uM)
            potassium_conc: float = 1.0,    # concentration of potassium (mM)
            magnesium_conc: float = 1.5,    # concentration of magnesium (mM)

            results_list = 10,              # maximum number of top results to show, best keep this at 10

            genome = 'Bos taurus'):         # genome to query for ePCR

        self.sequence = str(sequence.seq)
        self.bisulfite = bisulfite

        # strand is only important if bisulfite converted. If bisulfite converted and no strand specified, default to 'p' strand
        if bisulfite and strand is not None:
            self.strand = strand
        elif bisulfite and strand is None:
            self.strand = 'p'
        else:
            self.strand = None

        # everything else
        self.forward_start = forward_start if forward_start else 0
        self.forward_end = forward_end
        self.reverse_start = reverse_start
        self.reverse_end = reverse_end if reverse_end else len(self.sequence)
        self.PCR_max_length  = PCR_max_length

        self.primer_conc = primer_conc
        self.potassium_conc = potassium_conc
        self.magnesium_conc = magnesium_conc

        self.results_list = results_list

        self.genome = genome
        self.summary_conditions = (
            self.sequence,
            self.bisulfite,
            self.strand,

            self.forward_start,
            self.forward_end,
            self.reverse_start,
            self.reverse_end,

            self.PCR_max_length,
            self.primer_conc,
            self.potassium_conc,
            self.magnesium_conc,

            self.results_list,

            self.genome)

        self.primers_list = []

    def __repr__(self):
        return f'< {self.strand} (start-{self.forward_end}) -- ({self.reverse_start}-end) >'

    class BiSearch_primer:
        """
        BiSearch primers
        """
        ePCR_form_mapping = (
            'fp',		# forward primer sequence
            'rp',		# reverse primer sequence
            'bis',		# bisulfite status
            'db',		# genome database
            'fpcrlen'	# max ePCR length
            )

        def __init__(
                self,
                strand,
                score,
                f_seq,
                f_length,
                f_GC,
                f_Tm,
                r_seq,
                r_length,
                r_GC,
                r_Tm,
                amp_start,
                amp_end,
                amp_length,
                CpG_count
                ):

            self.strand = strand
            self.score = score
            self.f_seq = f_seq
            self.f_length = f_length
            self.f_GC = f_GC
            self.f_Tm = f_Tm
            self.r_seq = r_seq
            self.r_length = r_length
            self.r_GC = r_GC
            self.r_Tm = r_Tm
            self.amp_start = amp_start
            self.amp_end = amp_end
            self.amp_length = amp_length
            self.CpG_count = CpG_count

            self.primer_binds = []

        def __repr__(self) -> None:
            return f'< str: {self.strand}, f: {self.f_seq}, r: {self.r_seq}, score: {self.score}, length: {self.amp_length}, CpG count: {self.CpG_count}>'

        def do_BiSearch_ePCR(self, max_pcr_len=1000) -> None:
            """
            asdf
            """
            # define the ePCR conditions

            # initialize browser instance for BiSearch ePCR
            browser = mechanicalsoup.StatefulBrowser(soup_config={'features': 'lxml'})
            browser.open("http://bisearch.enzim.hu/?m=genompsearch")
            browser.select_form('form[action="?run"]')

            # deals with the degenerate primers
            for r in (('Y', 'T'),
                    ('R', 'A')):
                self.f_seq = self.f_seq.replace(*r)
            for r in (('Y', 'T'),
                    ('R', 'A')):
                self.r_seq = self.r_seq.replace(*r)

            list(map(browser.form.set, self.ePCR_form_mapping, (self.f_seq, self.r_seq, True, 'Bos taurus', max_pcr_len)))
            browser.submit_selected()

            #  find primer binds based on their position on the website
            # --------------------------------------------------

            # initialize primer binds
            self.primer_binds = [
                0,	# forward primer binds on the plus strand
                0,	# reverse primer binds on the plus strand
                0,	# forward primer binds on the minus strand
                0]	# reverse primer binds on the minus strand

            # define relative position of primer binds on strands
            pos_primer_matches = browser.get_current_page().find_all(text=re.compile('found based'))

            # for each in primer binds, find the values
            for i in range(len(self.primer_binds)):
                self.primer_binds[i] = pos_primer_matches[i].strip()[:-55] if len(pos_primer_matches[i].strip()[:-55]) > 0 else 0


            #  find  number of PCR products on the sense and antisense strand
            # --------------------------------------------------

            # define relative position of the PCR product number depending on how many PCR products are generated in the sense and antisense strand
            self.pos_PCR_prod_p = 1 if str(browser.get_current_page().find_all('h3')[3].find_next_siblings()[2]).strip() == '<br/>' else 2
            self.pos_PCR_prod_n = 1 if str(browser.get_current_page().find_all('h3')[6].find_next_siblings()[2]).strip() == '<br/>' else 2

            # number of PCR products is equal to the number of elements in those tables.
            self.PCR_products_p = len(browser.get_current_page().find_all('h3')[3].find_next_siblings()[self.pos_PCR_prod_p].find_all(text=re.compile('len')))
            self.PCR_products_n = len(browser.get_current_page().find_all('h3')[6].find_next_siblings()[self.pos_PCR_prod_n].find_all(text=re.compile('len')))


            #  find the preferred amplicon on either the sense or antisense strand
            # --------------------------------------------------

            # initialize  preferred amplicon's chromosome number and genomic position
            self.amp_chromosome_num = self.amp_genomic_start = self.amp_genomic_end = 'N/A'

            # define  relative web position of the preferred amplicon depending on the strand
            pos_amp = None
            if self.strand == 'p' and self.PCR_products_p > 0:
                pos_amp = 3
                pos_amp_str = self.pos_PCR_prod_p

            elif self.strand == 'm' and self.PCR_products_n > 0:
                pos_amp = 6
                pos_amp_str = self.pos_PCR_prod_n

            if pos_amp is not None:
                i_containing_amplicon = [i for i, s in enumerate(browser.get_current_page().find_all('h3')[pos_amp].find_next_siblings()[pos_amp_str].find_all(text=re.compile('len'))) if str(self.amp_length) in s][0]
                web_container = str(browser.get_current_page().find_all('h3')[pos_amp].find_next_siblings()[pos_amp_str].find_all('a')[i_containing_amplicon]).split('?')[1].split(';')

                self.amp_chromosome_num = web_container[0][4:-4]
                self.amp_genomic_start = web_container[2][6:-4]
                self.amp_genomic_end = web_container[3].split('"')[0][4:]

            self.num_degen_bases = sum(map(self.f_seq.count, ['R','Y'])) + sum(map(self.r_seq.count, ['R','Y']))
            self.num_repeats = sum(map(self.f_seq.count, ['AA','TT', 'CC', 'GG'])) + sum(map(self.r_seq.count, ['AA','TT', 'CC', 'GG']))

            return (
                self.strand,
                self.score,

                self.f_seq,
                self.f_length,
                self.f_GC,
                self.f_Tm,

                self.r_seq,
                self.r_length,
                self.r_GC,
                self.r_Tm,

                self.amp_start,
                self.amp_end,
                self.amp_length,
                self.CpG_count,

                self.primer_binds[0],
                self.primer_binds[1],
                self.PCR_products_p,
                self.primer_binds[2],
                self.primer_binds[3],
                self.PCR_products_n,

                self.amp_chromosome_num,
                self.amp_genomic_start,
                self.amp_genomic_end,

                self.num_degen_bases,
                self.num_repeats
            )

    def do_BiSearch_primer_design(self) -> list:
        """
        Performs BiSearch primer design and retrieves results.
        """

        # initialize browser instance for BiSearch ePCR, select form, fill in form according to class variables
        browser = mechanicalsoup.StatefulBrowser(soup_config={'features': 'lxml'})
        browser.open("http://bisearch.enzim.hu/?m=search")
        browser.select_form('form[action="?run"]')
        list(map(browser.form.set, self.BiSearch_form_mapping, self.summary_conditions, (True,)*len(self.summary_conditions)))
        browser.submit_selected()

        # translate HTML into something that pandas can read -- essentially gets rid of unnecessary form settings
        HTML_string = str(browser.get_current_page().find_all('table')[11])
        for r in (('<form action="?run" method="post">', ''),
                ('</form>', ''),
                ('<input name="prg" type="hidden" value="cgi/fpcr.cgi"/>', ''),
                ('<tr class="r1">', '<tr>'),
                ('<tr class="r2">', '<tr>')):
            HTML_string = HTML_string.replace(*r)
        BiSearch_Results = pd.read_html(HTML_string, skiprows=1)[0]

        # check if the results page is a legitimate table of primer results -- if there aren't 14 columns, there weren't any results
        if BiSearch_Results.shape[1] > 1:
            for row_i, _ in BiSearch_Results[::2].iterrows():
                x = BiSearch_Results.iloc

                # This is just some webscraping laziness
                score = x[row_i, 1]						# primer score, as scored by BiSearch
                strand = self.strand                    # strand
                f_seq = x[row_i, 2]                     # sequence of the foward primer
                f_length = len(x[row_i, 2])             # length of the foward primer
                f_GC = x[row_i, 5]                      # %GC of the foward primer
                f_Tm = x[row_i, 6]                      # T_m of the foward primer

                r_seq = x[row_i + 1, 2]                 # sequence of the reverse primer
                r_length = len(x[row_i + 1, 2])         # length of the reverse primer
                r_GC = x[row_i + 1, 5]                  # %GC of the reverse primer
                r_Tm = x[row_i + 1, 6]                  # T_m of the reverse primer

                amp_start = x[row_i, 3]                 # start of amplicon
                amp_end = (x[row_i + 1, 3])+2           # end of amplicon
                amp_length = (amp_end - amp_start)-1    # length of amplicon
                CpG_count = int(x[row_i, 8])            # number of CpG

                self.primers_list.append(self.BiSearch_primer(strand, score, f_seq, f_length, f_GC, f_Tm, r_seq, r_length, r_GC, r_Tm, amp_start, amp_end, amp_length, CpG_count))

        browser.close()
        return self.primers_list

# --------------------------------------------------
def main() -> None:
    """ Do the thing """

    args = get_args()

    fasta_file_path_arg = args.fasta_file_path.resolve()
    output_dir_arg, bis_arg, strand_arg, max_arg, sub_arg, tries_arg  = args.output_dir, args.bis, args.strand, args.max, args.sub, args.tries

    verbose_arg, seed_arg = args.verbose, args.seed
    random.seed(seed_arg)

    print_runtime('Started job')

    # check output_dir flag, if true: make an output dir
    if output_dir_arg:
        output_flag = 'BiSearch_primer_design'
        output_path_string = f'{output_flag}-output-{time.strftime("%Y_%m_%d_%H%M%S", time.localtime(time.time()))}'
        output_path = fasta_file_path_arg.joinpath(output_path_string)
        output_path.mkdir(parents=True, exist_ok=True)
    else:
        output_path = fasta_file_path_arg.parent

    # generate params
    generate_params_file(output_path, get_args())
    if verbose_arg > 0 : print_runtime('Created params.txt file')

    # generate and do
    list_of_fasta_files = []

    if fasta_file_path_arg.is_file():
        list_of_fasta_files = [fasta_file_path_arg]
    elif fasta_file_path_arg.is_dir():
        list_of_fasta_files = [fasta_file for fasta_file in fasta_file_path_arg.glob('*.fasta')]

    compiled_primers_total = {}
    for fasta_file in list_of_fasta_files:
        if verbose_arg > 0 : print_runtime(f'Generating BiSearch conditions for {fasta_file.name}')

        query_sequence = SeqIO.read(fasta_file, "fasta")
        fasta_subseqs = generate_subseqs(query_sequence, sub_arg, max_arg+200)

        subsequence_list = []
        for subseq_i in fasta_subseqs.keys():
            if tries_arg > len(fasta_subseqs[subseq_i].values()):
                if verbose_arg > 1 : print_runtime(f'!!! WARNING !!! : \nSplitting the FASTA into {tries_arg} tries in {len(fasta_subseqs.keys())} x {max_arg+200} bp subsequences isn\'t useful. The query region will overlap into the other regions anyway.')
                while tries_arg > len(fasta_subseqs[subseq_i].values()):
                    tries_arg -= 1
                if verbose_arg > 1 : print_runtime(f'Rounded down to {tries_arg} tries in {len(fasta_subseqs.keys())} regions to avoid this.')
            subsequence_list += random.choices(list(fasta_subseqs[subseq_i].values()), k=tries_arg)
        if verbose_arg > 1 : print_runtime(f'Generated {len(subsequence_list)} stepwise subsequence(s) to try')

        list_of_conditions = []
        for subsequence in subsequence_list:
            if verbose_arg > 1 : print_runtime('Sending conditions to BiSearch (could take up to 15 min...)')
            if bis_arg and not strand_arg:
                list_of_conditions.append(BiSearch_condition(
                        sequence = query_sequence,
                        reverse_start = subsequence['start'],
                        forward_end = subsequence['end'],
                        bisulfite = bis_arg,
                        strand = 'p',
                        ).do_BiSearch_primer_design()
                    )
                list_of_conditions.append(BiSearch_condition(
                        sequence = query_sequence,
                        reverse_start = subsequence['start'],
                        forward_end = subsequence['end'],
                        bisulfite = bis_arg,
                        strand = 'm',
                        ).do_BiSearch_primer_design()
                    )
            else:
                list_of_conditions.append(BiSearch_condition(
                        sequence = query_sequence,
                        reverse_start = subsequence['start'],
                        forward_end = subsequence['end'],
                        bisulfite = bis_arg,
                        strand = strand_arg,
                        ).do_BiSearch_primer_design()
                    )

            if verbose_arg > 1 : print_runtime('Received results from BiSearch')

        compiled_primers_for_fasta = []
        for primer_condition in list_of_conditions:
            for BiSearch_primer in primer_condition:
                if verbose_arg > 2 : print_runtime(f'Performed BiSearch ePCR with {BiSearch_primer}')
                compiled_primers_for_fasta.append(BiSearch_primer.do_BiSearch_ePCR())

        compiled_primers_total.update({fasta_file.stem: compiled_primers_for_fasta})

    for fasta in compiled_primers_total:
        generate_primers_csv(compiled_primers_total[fasta], fasta, output_path)
        print_runtime(f'Output .csv files containing primers for {fasta}.fasta')
    
    print_runtime(f'Finished job')

def generate_primers_csv(data, csv_name, output_path_arg) -> None:
    """
    Function outputs BiSearch-generated primers to .csv
    """

    pd.DataFrame(data,
            columns=[
                'strand',           # primer strand, 'p' = sense, 'm' = antisense
                'score',            # primer score, as scored by BiSearch

                'f_seq',            # sequence of forward primer
                'f_length',         # length of forward primer
                'f_GC',             # %GC of forward primer
                'f_Tm',             # T_m of forward primer

                'r_seq',            # sequence of reverse primer
                'r_length',         # length of reverse primer
                'r_GC',             # %GC of reverse primer
                'r_Tm',             # T_m of reverse primer

                'amp_start',        # start of amplicon
                'amp_end',          # end of amplicon
                'amp_size',         # size of amplicon
                'CpG_count',        # number of CpG

                'f_bind (+)',       # number of forward primer binds on (+) strand
                'r_bind (+)',       # number of reverse primer binds on (+) strand
                'PCR_products (+)', # number of PCR products on (+) strand

                'f_bind (-)',       # number of forward primer binds on (-) strand
                'r_bind (-)',       # number of reverse primer binds on (-) strand
                'PCR_products (-)', # number of PCR products on (-) strand

                'chr',              # chromosome number
                'start',            # start pos in chromosome
                'end',              # end pos in chromosome

                'degen_bases',      # number of degenerate bases in the primer sequences
                'repeats',          # number of repeats in the primer sequences
                ]).to_csv(output_path_arg.joinpath(f'{csv_name}.csv'))

def generate_subseqs(sequence, subseq_num, frame_length) -> dict:
    """
    Function will find subsequences of a given length within a fasta file.
    """
    subseq_len = int(len(sequence)/subseq_num)
    subseq_data = {}
    frame_length = subseq_len if frame_length > subseq_len else frame_length

    for subseq_count in range(subseq_num):
        subseq_list = {}

        subseq_start = subseq_count * subseq_len
        subseq_end = (subseq_count + 1) * subseq_len

        i = subseq_start
        while (i + frame_length) <= subseq_end:
            i += 1
            subseq_list.update({i: {
            'start': i,
            'end': i + frame_length,
            'CpG_count': sequence.seq[i:i+frame_length].count('CG')}})

        subseq_data.update({subseq_count: subseq_list})

    return subseq_data

def generate_params_file(output_path, args) -> None:
    """ Generates parameter files. """
    with open(output_path.joinpath(f'params-{time.strftime("%Y_%m_%d_%H%M%S", time.localtime(time.time()))}.txt'), 'w', encoding='UTF-8') as params_file:
        output_vars = [
            f'input={args.fasta_file_path.resolve()}',
            f'output_dir={args.output_dir}',
            f'bis={args.bis}',
            f'strand={args.strand}',
            f'max_PCR_len={args.max}',
            f'n_subseq={args.sub}',
            f'n_tries={args.tries}',
            f'verbose={args.verbose}',
            f'seed={args.seed}']
        params_file.write('\n'.join(output_vars))

def print_runtime(action) -> None:
    """ Print the time and some defined action. """
    print(f'[{time.strftime("%Y-%m-%d %H:%M:%S", time.localtime())}] {action}')

# --------------------------------------------------
if __name__ == '__main__':
    main()