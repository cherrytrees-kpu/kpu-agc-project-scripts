from Bio import SeqIO, Seq
import pathlib
import argparse
import tempfile
import io
import subprocess
import pandas
import csv
import multiprocessing

class multi_blast(): 
    def __init__(self, query_seq): 
        self.query = query_seq
    def _blast(seq, blastdb, blastdb_len): 
        #Create a temporary file containing the sequence to pass as an input to BLAST
        #Temporary file will be deleted after the program finishes running
        #Sequence is encoded to utf-8, and we move to start of temporary file
        fasta = tempfile.NamedTemporaryFile(delete=True)
        fasta.write(f">probe\n{str(seq)}".encode())
        fasta.seek(0)
        #blastn command
        args = [
            "blastn",
            "-task",
            "blastn-short",
            "-db",
            blastdb,
            "-num_alignments",
            str(blastdb_len),
            "-outfmt",
            "10 qacc sacc ssciname pident qlen length mismatch gapopen qstart qend sstart send evalue bitscore",
            "-query",
            fasta.name,
        ]
        #Run the program and decode output
        #Turn io stream into string
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
    def run(self): 
        #Split the job
        args = []
        for i in self.query: 
            args.append((i, '/data/db/SSU-nematodes.fasta', '17095'))
        with multiprocessing.Pool() as p: 
            blast_results = p.starmap(_blast, args)
            for result in blast_results: 
                print(result)

def main(): 
    seq = ['atcgatgcatgatc', 'atgatatagatcggt', 'atgatcggatagcagat', 'atagcatagagtaca']
    bob = multi_blast(seq)
    bob.run()

if __name__ == '__main__': 
    main()