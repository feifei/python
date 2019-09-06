'''
    Analyze the ASH level of different genomes
    To have a consistent estimation to compare
    BWA -> pileup -> snps -> calculation
    First steps in sh, easier to fix manually
'''
import re, os
import argparse
from collections import defaultdict
from Bio.SeqUtils import GC
from Bio import SeqIO
from Bio.Seq import Seq


def clean_seq(seq):
    s = re.sub("[+-]\d+[A-Za-z]+", "", seq)
    s = re.sub("\^.|\$|\,|\.", "", s)
    return s


def get_count(seq):
    d = defaultdict(int)
    if not seq:
        return None
    d["A"] = seq.count("A") + seq.count("a")
    d["C"] = seq.count("C") + seq.count("c")
    d["G"] = seq.count("G") + seq.count("g")
    d["T"] = seq.count("T") + seq.count("t")
    
    return d

# Make it part of input 0.1 and 20
def get_alt_bases(count_d, cov):
    s = ""
    for k, v in count_d.iteritems():
        p = v/float(cov)
        if p >= 0.1 and cov >= 20:
            s += "%s:%d:%.2f;" %(k, v, p)
    return s


parser = argparse.ArgumentParser(description='Estimate ASH')
parser.add_argument('reference', type=argparse.FileType('r'))
parser.add_argument('pileup_file')
parser.add_argument('--cov', dest='cov_cutoff', default=20, type=int,
                    help='Site coverage cutoff')
parser.add_argument('--percent', dest='percent', default=0.1, type=float,
                    help='Percent of alternative base cutoff')

args = parser.parse_args()
pileup_file = args.pileup_file
reference = args.reference
basename, extension = os.path.splitext(args.pileup_file)
snps_outfile = basename + ".snps.tab"

nr_snps = 0
with open(snps_outfile, 'w') as outh, open(pileup_file) as pileup_inh:
    for line in pileup_inh:
        arr = line.split("\t")
        if len(arr) != 6:
            print line
            continue
        scfid, pos, ref_base, cov, seq, qual = arr
        c_seq = clean_seq(seq)
        count_d = get_count(c_seq)
        if not count_d:
            continue
        alt_bases = get_alt_bases(count_d, int(cov))
        if not alt_bases:
            continue
        
        nr_snps += 1
        print >>outh, "\t".join([scfid, pos, ref_base, cov, alt_bases])
    

seq_d = defaultdict(Seq)
tot_len = 0
tot_len_without_N = 0
for record in SeqIO.parse(reference, 'fasta'):
    seq_d[record.id] = record.seq
    tot_len += len(record.seq)
    tot_len_without_N += len(re.sub("N+", "", str(record.seq).upper()))

mean_GC = GC("".join(map(str, seq_d.values())))
ash = nr_snps / float(tot_len)
ash_without_N = nr_snps / float(tot_len_without_N)

print "Genome size: %d; Genome size without N %d; in %d contigs" %(tot_len, tot_len_without_N, len(seq_d))
print "ASH %.3f%%; ASH (genome size without N) %.3f%%; GC %.2f%%" %(ash*100, ash_without_N*100, mean_GC)

