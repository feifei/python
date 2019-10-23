#!/usr/bin/env python
'''
    Extend blast hit to ORF
'''

from Bio import SeqIO
from translate_l import *
from frame import *


def get_write_seq(fasta_file, blast, out):
    seq_dict = {rec.id : rec.seq for rec in SeqIO.parse(fasta_file, "fasta")}
    aa_out = out + ".faa"
    na_out = out + ".fa"
    i = 1
    with open(blast, 'r') as inh, open(aa_out, 'w') as aa_fh, open(na_out, 'w') as na_fh:
        for line in inh:
            line = line.strip()
            if not line or len(line) == 0:
                continue
        
            arr = line.split("\t")
            scfid = arr[0]
            start, end, frame, ref_length = map(int, arr[1:-1])
            evalue = float(arr[-1])
            if start > end:
                start, end = end, start
            if (end - start) / 3 < ref_length * 0.9:
                continue
            seq = seq_dict[scfid]                
            start_f = normal_pos_to_frame(seq, frame, start)
            end_f = normal_pos_to_frame(seq, frame, end)
            seq_f = seq.frame(frame)
            x = longest_orf_nearby(seq_f, ref_length, start_f, end_f, table=6, stop_symbol="*", start_symbol="M", to_stop=False)
            if x:
                na_seq, aa_seq, start, end = x
                print >>na_fh, ">%d" %i
                print >>na_fh, na_seq
                print >>aa_fh, ">%d" %i
                print >>aa_fh, aa_seq
                i += 1




# vortens
fasta_file = "/Users/feifei/sequenceserver_db/pb_114_polished_assembly.fasta"
blast = "/Users/feifei/Projects/Others/Asgeir/vortens.blast.out"
out = "/Users/feifei/Projects/Others/Asgeir/CHC.vortens"
get_write_seq(fasta_file, blast, out)


# barkhanus pb_116
fasta_file = "/Users/feifei/sequenceserver_db/pb_116_polished_assembly.fasta"
blast = "/Users/feifei/Projects/Others/Asgeir/barkhanus.blast.out"
out = "/Users/feifei/Projects/Others/Asgeir/CHC.barkhanus"
get_write_seq(fasta_file, blast, out)

# hexamita pb_188
fasta_file = "/Users/feifei/sequenceserver_db/pb_188_polished_assembly.fasta"
blast = "/Users/feifei/Projects/Others/Asgeir/hexamita.blast.out"
out = "/Users/feifei/Projects/Others/Asgeir/CHC.hexamita"
get_write_seq(fasta_file, blast, out)



