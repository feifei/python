'''
    Extend blast hit to ORF
'''

from Bio import SeqIO
from translate_l import *
from frame import *


fasta_file = "/Users/feifei/sequenceserver_db/pb_114_polished_assembly.fasta"

# matches info:
matches = \
"""
unitig_253|quiver   112038  111223  -1  427

unitig_4094|quiver  61686   62438   3   511

unitig_10585|quiver 37462   38382   1   467

unitig_10585|quiver 34666   35745   1   511

unitig_10585|quiver 36034   37266   1   404

unitig_8882|quiver   16592  15180  -2  501

unitig_8882|quiver  16568   15402   -2  468
"""

seq_dict = {rec.id : rec.seq for rec in SeqIO.parse(fasta_file, "fasta")}

for line in matches.split("\n"):
    if not line or len(line) == 0:
        continue
        
    arr = line.split()
    scfid = arr[0]
    print arr
    start, end, frame, ref_length = map(int, arr[1:])
    seq = seq_dict[scfid]
    if start > end:
        start, end = end, start
    start_f = normal_pos_to_frame(seq, frame, start)
    end_f = normal_pos_to_frame(seq, frame, end)
    seq_f = seq.frame(frame)
    na_seq, aa_seq, start, end = longest_orf_nearby(seq_f, ref_length, start_f, end_f, table=6, stop_symbol="*", start_symbol="M", to_stop=False)
    
    print scfid, aa_seq
