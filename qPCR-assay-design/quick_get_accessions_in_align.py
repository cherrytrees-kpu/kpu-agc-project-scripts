import sys
import pathlib
from Bio import AlignIO

path = pathlib.Path(sys.argv[1])

alignment = AlignIO.read(path, 'fasta')

list_id = []
for seq in alignment: 
    list_id.append(seq.id.split('.')[0])

output_file = open('list_accessions.txt', 'w')
for id in list_id: 
    output_file.write(id+'\n')
output_file.close()