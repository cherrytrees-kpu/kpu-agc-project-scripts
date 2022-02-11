import vcf
import pathlib
import argparse
import collections

class generateVCF: 
    def __init__(self, vcf_path, output_path, id): 
        self.vcf_directory = vcf_path
        self.output_path = output_path
        self.id = id
        self.records = []
        self.vcf_reader_template = None
    def retrieve_records(self): 
        def check_record(record): 
            if(
                len(record.FILTER) > 0
                or record.num_het > 0
            ):
                return False
            else: 
                return True
        for vcf_path in self.vcf_directory.glob('*.vcf'): 
            print(str(vcf_path))
            vcf_file = open(vcf_path, 'r')
            vcf_reader = vcf.Reader(vcf_file)
            for record in vcf_reader: 
                if check_record(record) is True: 
                    self.records.append(record)
    def create_het_records(self): 
        CallData = collections.namedtuple('CallData', ['GT', 'GQ'])
        for record in self.records:
            #Add REF to ALT list
            record.samples[0].data = CallData('0/1', 0)
    def output_vcf(self): 
        #Generate template
        template_path = next(self.vcf_directory.glob('*.vcf'))
        template = vcf.Reader(open(template_path, 'r'))
        #Create output
        output_file_path = self.output_path.joinpath(f'{self.id}_genotypes.vcf')
        output_file = open(output_file_path, 'w')
        vcf_writer = vcf.Writer(output_file, template)
        for record in self.records:
            print(record) 
            vcf_writer.write_record(record)
        output_file.close()

def parse_args(): 
    parser = argparse.ArgumentParser("Program to prepare a VCF file for NGS read simulation")
    parser.add_argument(
        'vcf_path', 
        type=pathlib.Path, 
        action='store', 
        help='Path to directory containing VCF files'
    )
    parser.add_argument(
        'id', 
        type=str,
        action='store', 
        help='ID to name the output VCF file'
    )
    parser.add_argument(
        '-o',
        '--output',
        dest='output_path',
        type=pathlib.Path,
        action='store', 
        default=None,
        help='Output path. Default: output to where the VCF path is.'
    )
    args = parser.parse_args()
    if args.output_path is None: 
        args.output_path = args.vcf_path
    return (args.vcf_path, args.output_path, args.id)

def main(): 
    vcf_path, output_path, id = parse_args()
    prep_vcf = generateVCF(vcf_path, output_path, id)
    prep_vcf.retrieve_records()
    #CURRENTLY, SIMUG JUST IGNORES HETEROZYGOUS RECORDS... I WANT TO VOMIT
    #prep_vcf.create_het_records()
    prep_vcf.output_vcf()

if __name__ == '__main__': 
    main()