#!/usr/bin/python

import sys
import subprocess
from optparse import OptionParser
import re


def main():
    parser = OptionParser()
    parser.add_option("-i", "--infile", dest="infile",
                      help="fasta format of NT/AA alignment sequences", metavar="FILE")
    parser.add_option("-d", "--work_dir", dest="work_dir",
                      default=".",
                      help="working directory")
    
    
    (options, args) = parser.parse_args()
    
    infile = options.infile
    work_dir = options.work_dir
    if infile == None:
        print "use -h or --help to see the help"
        sys.exit(1)
        
    if work_dir == None:
        run_geneconv(infile)
    else:
        run_geneconv(infile, work_dir)
    

def run_geneconv(infile, cfgfile, work_dir='.'):
    cl = GeneconvCommandline(infile, cfgfile)
    child = subprocess.Popen(str(cl), 
                     stdout=subprocess.PIPE,
                     stderr=subprocess.PIPE,
                     cwd=work_dir,
                     shell=True)
    
    # Close pipes to prevent them from filling up, since we're not reading from them.
    child.stdout.close()
    child.stderr.close()
    
    # wait for the subprocess to finish first
    child.wait()
    #assert child.returncode == 0, "geneconv failed"

# Check if there's GI (internal recom) or GO (outer-seq fragment) 
def parseFragsFiles(frags_file):
    frags_handle = open(frags_file, 'r')
    flag = False
    for line in frags_handle:
        if re.match('GI', line):
            flag = True
        elif re.match('GO', line):
            flag = True
    return flag


# Return the GI and GO results in table format, every row represents a GI 
def fragToTable(frags_file):

    if parseFragsFiles(frags_file):
        row = []
        with open(frags_file, 'r') as frags_h:
            for line in frags_h:
                if re.match('[GI|GO]', line):
                    arr = line.split()
                    tag = arr[0]
                    geneids = arr[1].split(";")
                    geneid1 = geneids[0]
                    if len(geneids) != 2:
                        geneid2 = None
                    else:
                        geneid2 = geneids[1]                        
                    sim_pval = float(arr[2])
                    bc_ka_pval = float(arr[3])
                    begin = int(arr[4])
                    end = int(arr[5])
                    length = int(arr[6])
                    poly_num = int(arr[7])
                    diff_num = int(arr[8])
                    tot_diff = int(arr[9])
                    mism_pen = int(arr[10])
                    row.append([tag, geneid1, geneid2, sim_pval, bc_ka_pval, begin, end, length, poly_num, diff_num, tot_diff, mism_pen])
        
        return row
                    
    else:
        return None


class GeneconvCommandline:
    """ Commandline for geneconv.
    """
    def __init__(self, infile, cfgfile):
        self.infile = infile
        self.cfgfile = cfgfile
        # self.parameters = []
        
    def __str__(self):
        return "geneconv %s %s" % (self.infile, self.cfgfile)
    
    
if __name__ == "__main__":
    main()
