#!/usr/bin/env python

import os
import subprocess
from optparse import OptionParser
from Bio import SeqIO


def main():

    parser = OptionParser("usage: %prog fasta_file num_threads ")

    (options, args) = parser.parse_args()

    if len(args) != 2:
        parser.error("incorrect number of arguments")

    fasta_file = args[0]
    num_threads = int(args[1])
    split_fasta_file_by_num_threads(fasta_file, num_threads)
    
    
def split_fasta_file_by_num_threads(fasta_file, num_threads):
    file_name, file_ext = os.path.splitext(fasta_file)
    tot_num_records = int(subprocess.check_output(["grep", "-c", "^>", fasta_file]))
    num_records = tot_num_records / num_threads
    
    count = 0
    records = []
    file_num = 0
    with open(fasta_file, 'r') as inh:
        for record in SeqIO.parse(inh, 'fasta'):
            count += 1
            records.append(record)
            if count % num_records == 0:
                file_num += 1
                outfile = "%s.%d%s" %(file_name, file_num, file_ext)
                with open(outfile, 'w') as outh:
                    SeqIO.write(records, outh, 'fasta')
                records = []
                count = 0

    if records:
        outfile = "%s.%d%s" %(file_name, file_num, file_ext)
        with open(outfile, 'a') as outh:
            SeqIO.write(records, outh, 'fasta')

def remove_splitted_fasta_files(fasta_file):
    file_name, file_ext = os.path.splitext(fasta_file)
    subprocess.call(["rm", "%s.*%s" %(file_name, file_ext)])        


if __name__ == "__main__":
    main()
