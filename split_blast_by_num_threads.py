''' 
    Splitting blast job according to the specified number of threads 
    Only blastp for now
'''
import threading

from blast import *
from split_fasta_file_by_num_threads import *


class blastpThread (threading.Thread):
    def __init__(self, infile, db, evalue, outfmt, outfile):
        threading.Thread.__init__(self)
        self.infile = infile
        self.db = db
        self.evalue = evalue
        self.outfmt = outfmt
        self.outfile = outfile
    def run(self):
        run_blastp(self.infile, self.db, args = {}, evalue = self.evalue, outfmt = self.outfmt, out= self.outfile)


def main():
    parser = OptionParser("usage: %prog fasta_file num_threads db outfile")
    parser.add_option("-e", "--evalue", dest="evalue", default=10, type=float,
                      help="Spcify evalue cutoff")
    parser.add_option("-f", "--outfmt", dest="outfmt", default = 0, type=int,
                      help="Spcify blast output format")
    
    
    (options, args) = parser.parse_args()

    if len(args) != 4:
        parser.error("incorrect number of arguments")

    fasta_file = args[0]
    infile_name, infile_ext = os.path.splitext(fasta_file)
    num_threads = int(args[1])
    db = args[2]
    blast_outfile = args[3]
    outfile_name, outfile_ext = os.path.splitext(blast_outfile)
    evalue = options.evalue
    outfmt = options.outfmt
    
    split_fasta_file_by_num_threads(fasta_file, num_threads)
    
    # Run blast mutltithreading
    threads = []
    for i in range(num_threads):
        infile = "%s.%d%s" %(infile_name, i + 1, infile_ext)
        outfile = "%s.%d%s" %(outfile_name, i + 1, outfile_ext)
        print infile, outfile
        thread = blastpThread(infile, db, evalue, outfmt, outfile)
        thread.start()
        threads.append(thread)
    
    # Wait for all threads to complete
    for t in threads:
        t.join()

    # Merge blast outputs
    with open(blast_outfile, 'w') as outh:
        for i in range(num_threads):
            infile = "%s.%d%s" %(infile_name, i + 1, infile_ext)
            outfile = "%s.%d%s" %(outfile_name, i + 1, outfile_ext)
            with open(outfile, 'r') as inh:
                outh.write(inh.read())
    
            
            # Deleting tmp files
            subprocess.call(['rm', outfile, infile])
    


if __name__ == "__main__":
    main()


