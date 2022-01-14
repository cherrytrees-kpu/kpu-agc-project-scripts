#!/usr/bin/env python3
"""
Author : Erick Samera
Date   : 2022-01-08
Purpose: Calculates the Tm of primers and estimates an appropriate annealing temp. for different polymerases.
"""

import argparse
from typing import NamedTuple

import pathlib
import math

import pandas as pd
import time

class Args(NamedTuple):
    """ Command-line arguments """
    file: pathlib.Path
    conc: float
    pol: str

# --------------------------------------------------
def get_args() -> Args:
    """ Get command-line arguments """

    parser = argparse.ArgumentParser(
        description='Batch ThermoFisher Tm Calculator',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument('file',
                        metavar='file',
                        type=pathlib.Path,
                        help='Path of Thermo-formatted .txt file containing primer list')

    parser.add_argument('-c',
                        '--conc',
                        help="Primer concentration (uM)",
                        metavar='conc',
                        type=float,
                        default=0.5)

    parser.add_argument('-p',
                        '--pol',
                        help="Specify polymerase to use {SuperFi, Phusion, DreamTaq}",
                        metavar='pol',
                        type=str,
                        choices=['SuperFi', 'Phusion', 'DreamTaq'],
                        default='Phusion')

    args = parser.parse_args()

    if not args.file.is_file and args.file.suffix == '.txt':
        parser.error(f'Specify a .txt file')

    return Args(args.file, args.conc, args.pol)

# --------------------------------------------------
class primer:

    def __init__(
            self,
            name,
            seq,
            primer_conc,
            ):
        self.name = name
        self.seq = seq.upper()
        self.primer_conc = primer_conc * 1e-6
        self.salt_conc = 50 * 1e-3
        return None
    
    def __repr__(self):
        return f"<{self.name} '{self.seq}'>"

    def nucleotide_count(self):
        """
        Function will get count of each nucleotide in primer.
        """

        nucleotides = {
            'A': self.seq.count('A'),
            'T': self.seq.count('T'),
            'C': self.seq.count('C'),
            'G': self.seq.count('G'),
            }

        return nucleotides
    
    def GC_content(self):
        """
        Function will get GC content of primer.
        """

        G_count = self.seq.count('G')
        C_count = self.seq.count('C')
        return (G_count + C_count)/len(self.seq)*100

    def molecular_weight(self):
        """
        Function will get molecular weight of primer.
        """

        molecular_weight_DNA_data = {
            'A': 313.21,
            'T': 304.2,
            'C': 289.18,
            'G': 329.21,
            'I': 314
            }

        molecular_weight = 0

        for base in self.seq:
            molecular_weight += molecular_weight_DNA_data[base]

        return molecular_weight-61.96
    
    def sumObjectSymbolNNValues(self, values, inits):
        """
        Get the dS and dH from primer sequence using different methods.
        """

        total_dH = inits['dH_ini']
        total_dS = inits['dS_ini']

        for i, s in enumerate(self.seq):
            if i < len(self.seq)-1:
                total_dH += values[s][self.seq[i+1]]['dH']
                total_dS += values[s][self.seq[i+1]]['dS']
        
        return {'dS': total_dS, 'dH': total_dH}

    def Tm_simple(self):
        """
        Get the Tm using simple calculation.
        """
        
        Tm_simple_params_short = {
            'A': 2,
            'T': 2,
            'C': 4,
            'G': 4
            }

        Tm_simple_params_long = {
            'A': 64.9,
            'T': 64.9,
            'G': 64.9 + 41,
            'C': 64.9 + 41
            }
        
        if len(self.seq) < 14:
            return sum([Tm_simple_params_short[nucleotide] * self.nucleotide_count()[nucleotide] for nucleotide in ['A', 'T', 'C', 'G']])
        
        else:
            return (sum([Tm_simple_params_long[nucleotide] * self.nucleotide_count()[nucleotide] for nucleotide in ['A', 'T', 'C', 'G']])- 41 *16.4)/len(self.seq)

    def Tm_All97(self):
        """
        Get the Tm using Allawi-SantaLucia (1997) method.
        """
        
        Tm_All97_data = {
            'A':{
                'A':{'dS':-22.2,'dH':-7900.0},
                'C':{'dS':-22.4,'dH':-8400.0},
                'T':{'dS':-20.4,'dH':-7200.0},
                'G':{'dS':-21.0,'dH':-7800.0}},
            'C':{
                'A':{'dS':-22.7,'dH':-8500.0},
                'C':{'dS':-19.9,'dH':-8000.0},
                'T':{'dS':-21.0,'dH':-7800.0},
                'G':{'dS':-27.2,'dH':-10600.0}},
            'T':{
                'A':{'dS':-21.3,'dH':-7200.0},
                'C':{'dS':-22.2,'dH':-8200.0},
                'T':{'dS':-22.2,'dH':-7900.0},
                'G':{'dS':-22.7,'dH':-8500.0}},
            'G':{
                'A':{'dS':-22.2,'dH':-8200.0},
                'C':{'dS':-24.4,'dH':-9800.0},
                'T':{'dS':-22.4,'dH':-8400.0},
                'G':{'dS':-19.9,'dH':-8000.0}}}

        dS_ini = 0
        dH_ini = 0

        if self.seq.endswith(('A', 'T')):
            dS_ini += 4.1
            dH_ini += 2300
        if self.seq.endswith(('C', 'G')):
            dS_ini -= 2.8
            dH_ini += 100
        if self.seq.startswith(('A', 'T')):
            dS_ini += 4.1
            dH_ini += 2300
        if self.seq.startswith(('C', 'G')):
            dS_ini -= 2.8
            dH_ini += 100

        thermodyndata = self.sumObjectSymbolNNValues(Tm_All97_data, {'dH_ini': dH_ini, 'dS_ini': dS_ini})
        res = (thermodyndata['dH'] / (1.9872 * math.log(self.primer_conc / 4.0) + thermodyndata['dS'])) + (16.6 * math.log(0.215273974689348) / math.log(10)) - 273.15

        res_adj = (res+3)*0.9376798568+4.5185404499

        if res_adj < 0:
            res_adj = 0
        elif res_adj > 95:
            res_adj =  95

        return res_adj

    def Tm_Breslauer(self):
        """
        Get the Tm using Breslauer method.
        """

        Tm_Breslauer_data = {
            'A':{
                'A':{'dH':-9100,'dS':-24},
                'T':{'dH':-8600,'dS':-23.9},
                'G':{'dH':-7800,'dS':-20.8},
                'C':{'dH':-6500,'dS':-17.3}},
            'T':{
                'A':{'dH':-6000,'dS':-16.9},
                'T':{'dH':-9100,'dS':-24},
                'G':{'dH':-5800,'dS':-12.9},
                'C':{'dH':-5600,'dS':-13.5}},
            'G':{
                'A':{'dH':-5600,'dS':-13.5},
                'T':{'dH':-6500,'dS':-17.3},
                'G':{'dH':-11000,'dS':-26.6},
                'C':{'dH':-11100,'dS':-26.7}},
            'C':{
                'A':{'dH':-5800,'dS':-12.9},
                'T':{'dH':-7800,'dS':-20.8},
                'G':{'dH':-11900,'dS':-27.8},
                'C':{'dH':-11000,'dS':-26.6}}}
        
        thermodyndata = self.sumObjectSymbolNNValues(Tm_Breslauer_data, {'dH_ini': -3400, 'dS_ini': -12.4})
        res = (thermodyndata['dH'] / (1.9872 * math.log(self.primer_conc / 160) + thermodyndata['dS'])) + (16.6 * math.log(self.salt_conc) / math.log(10)) - 273.15
        
        if res < 0:
            res = 0
        elif res > 95:
            res = 95
        return res

    def Tm_Taq(self):
        """
        Get the Tm using maybe the SantaLucia method (1996)?
        """

        Tm_San96_data = {
            'A':{
                'A':{'dS':-23.6,'dH':-8400.0},
                'C':{'dS':-23.0,'dH':-8600.0},
                'T':{'dS':-18.8,'dH':-6500.0},
                'G':{'dS':-16.1,'dH':-6100.0}},
            'C':{
                'A':{'dS':-19.3,'dH':-7400.0},
                'C':{'dS':-15.6,'dH':-6700.0},
                'T':{'dS':-16.1,'dH':-6100.0},
                'G':{'dS':-25.5,'dH':-10100.0}},
            'T':{
                'A':{'dS':-18.5,'dH':-6300.0},
                'C':{'dS':-20.3,'dH':-7700.0},
                'T':{'dS':-23.6,'dH':-8400.0},
                'G':{'dS':-19.3,'dH':-7400.0}},
            'G':{
                'A':{'dS':-20.3,'dH':-7700.0},
                'C':{'dS':-28.4,'dH':-11100.0},
                'T':{'dS':-23.0,'dH':-8600.0},
                'G':{'dS':-15.6,'dH':-6700.0}}}

        thermodyndata = self.sumObjectSymbolNNValues(Tm_San96_data, {'dH_ini': 0.0, 'dS_ini': -0.0})
        NaEquiv = 0.15527397
        nucleotide_F_term  = -15.894952
        etropyCorrection =  0.368 * (len(self.seq) - 1.0) * math.log(NaEquiv); 
        thermodyndata['dS'] += + etropyCorrection

        res = thermodyndata['dH'] / (thermodyndata['dS'] + 1.9872 * nucleotide_F_term) - 273.15
        
        if res < 0:
            res=0
        elif res > 95:
            res = 95
        
        res_ajusted = res * 0.9085395477132917 - 3.707388372789194 if (len(self.seq) < 21) else (res + 3) * 0.9085395477132917 - 3.707388372789194

        return res_ajusted

# --------------------------------------------------
def main() -> None:

    start_main = time.time()
    
    # define args
    args = get_args()
    file_arg = args.file
    conc_arg = args.conc
    pol_arg = args.pol

    primers_list = []

    # open primer file, get values, and run through calculate_Ta()
    with open(file_arg, 'r') as primer_file:
        for line in primer_file:

            # separate out lines by ';' and get the name, sequence values
            primer_1_val = line.split(';')[0].strip()
            primer_1 = primer(primer_1_val.split(' ')[0], primer_1_val.split(' ')[1], conc_arg)

            primer_2_val = line.split(';')[1].strip()
            primer_2 = primer(primer_2_val.split(' ')[0], primer_2_val.split(' ')[1], conc_arg)

            Ta, note  = calculate_Ta(primer_1, primer_2, pol_arg)

            primer_info = [
                primer_1.name,
                primer_1.seq, 
                round(calculate_Tm(primer_1, pol_arg), 1), 
                
                primer_2.name, 
                primer_2.seq, 
                round(calculate_Tm(primer_2, pol_arg), 1),

                round(Ta, 1),
                note
                ]

            primers_list.append(primer_info)
    
    pd.DataFrame(primers_list,
        columns=['Primer 1', 'Primer 1 Sequence', 'Tm', 'Primer 2', 'Primer 2 Sequence', 'Tm', 'Ta', 'notes']).to_csv(
            file_arg.parent.joinpath(file_arg.stem + '.csv'), 
            index=None)
    
    print_runtime(round(time.time() - start_main, 3), f'get calculate Ta information in {file_arg}')

def calculate_Tm(primer, pol_arg):
    """
    Function will calculate melting temperature (Tm) based on polymerase.
    """

    if pol_arg in ['SuperFi', 'Phusion']:
        Tm = primer.Tm_All97()

    elif pol_arg in ['DreamTaq']:
        Tm = primer.Tm_Taq()
    
    return Tm
            
def calculate_Ta(primer_1, primer_2, pol_arg):
    """
    Function will calculate annealing temperature (Ta) based on polymerase.
    """

    if pol_arg in ['SuperFi', 'Phusion']:
        Ta = min(primer_1.Tm_All97(), primer_2.Tm_All97()) if ((min(len(primer_1.seq), len(primer_2.seq)) < 21) and (min(primer_1.Tm_All97(), primer_2.Tm_All97()) < 72)) else min(min(primer_1.Tm_All97(), primer_2.Tm_All97()), 72)
    
    elif pol_arg in ['DreamTaq']:
        Ta = min((max(primer_1.Tm_Taq(), primer_2.Tm_Taq())), 72) if (min(primer_1.Tm_Taq(), primer_2.Tm_Taq()) < 72) else 72

    # notes to output in case of errors, maybe?
    # --------------------------------------------------
    note = ''

    if abs(calculate_Tm(primer_1, pol_arg) - calculate_Tm(primer_2, pol_arg)) >= 5:
        note += "Tm difference of more than 5C or greater is not recommended. "
    if (Ta < 45):
        note += "Annealing temperature lower than 45C is not recommended. "
    if (len(primer_1.seq) < 7) or (len(primer_2.seq) < 7):
        note += "Both primers need to be longer than 7 nt. "
    if ((calculate_Tm(primer_1, pol_arg) > 69) and (calculate_Tm(primer_1, pol_arg) < 72)) and ((calculate_Tm(primer_2, pol_arg) > 69) and (calculate_Tm(primer_2, pol_arg) < 72)):
        note += "A 2-step protocol (combined annealing/extension) is recommended when primer Tm values are higher than 69C, using 72C for annealing step. "
    if (Ta >= 72):
        note += "Annealing temperature should not exceed 72C. "
    return Ta, note

def print_runtime(runtime, action):
    print(f'--- Took {runtime} s to {action}. ---')

# --------------------------------------------------
if __name__ == '__main__':
    main()
