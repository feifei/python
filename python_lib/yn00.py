import sys
import subprocess
from optparse import OptionParser


def main():
	parser = OptionParser()
	parser.add_option("-i", "--infile", dest="infile",
	                  help="fasta format of DNA sequences aligned", metavar="FILE")
	parser.add_option("-o", "--outfile_aln", dest="outfile", default = "nuc_subs.yn",
	                  help="result file from yn00 program, default is nuc_subs.yn", metavar="FILE")
	parser.add_option("-d", "--directory", dest="dir", default = ".",
	                  help="working directory, default is the current directory")
	
	
	(options, args) = parser.parse_args()

	infile = options.infile
	outfile = options.outfile
	work_dir = options.dir
	if infile == None or outfile == None:
		print "use -h or --help to see the help"
		sys.exit(1)

	fixed_infile = fix_file_file(infile, infile+".fixed")
	run_yn00(fixed_infile, outfile, work_dir)
	
	
def find_synonymous(yn_file="nuc_subs.yn", work_dir="."):
		
	ds_values = []
	dn_values = []
	seq_num1s = []
	seq_num2s = []

	output_h = open(work_dir + "/" + yn_file)
	for line in output_h.xreadlines():
		if line.find("+-") >= 0 and line.find("dS") == -1:
			ds_value, dn_value, seq_num1, seq_num2 = extract_values(line)
			ds_values.append(ds_value)
			dn_values.append(dn_value)
			seq_num1s.append(seq_num1)
			seq_num2s.append(seq_num2)

	
	if ds_values == [] or dn_values == []:
		print "yn00 didn't work: \n%s" % output_h.read()
	
	return ds_values, dn_values, seq_num1s, seq_num2s



def run_yn00(infile, outfile="nuc_subs.yn", ctl_file = "yn-input.ctl", work_dir="."):	
	"""Run yn00 to find the synonymous subsitution rate for the alignment.
	"""
	# create the .ctl file
	ctl_h = open(work_dir + "/" + ctl_file, "w")
	ctl_h.write("seqfile = %s\noutfile = %s\nverbose = 0\n" % 
	            (infile, outfile))
	ctl_h.write("icode = 0\nweighting = 0\ncommonf3x4 = 0\n")
	ctl_h.close()
	
	cl = YnCommandline(ctl_file)
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
	assert child.returncode == 0, "yn00 failed"


# ouput records in yn00's format, sort of phylip format
def fix_file_records(records, outfile):
	output_h = open(outfile, "w")

	output_h.write("\t%d %d\n" % (len(records), len(records[0].seq)))
	for record in records:
		output_h.write(">"+record.id+"\n")
		output_h.write(str(record.seq)+"\n")
	output_h.close()
	return outfile

def fix_file_file(infile, outfile):
	
	input_h = open(infile)
	records = list(SeqIO.parse(fh, "fasta"))
	output_h = open(outfile, "w")

	output_h.write("\t%d %d\n" % (len(records), len(records[0].seq)))
	for record in records:
		output_h.write(">"+record.id+"\n")
		output_h.write(str(record.seq)+"\n")
	output_h.close()
	
	return outfile

	
def extract_values(line):
	"""Extract subsitution values from a line of text with the information of sequence number.

	This is just a friendly function to grab a float value for Ks and Kn 
	values from the junk I get from the last line of the yn00 file.

	Line:
	2    1    52.7   193.3   2.0452  0.8979  0.0193 0.0573 +- 0.0177 
	2.9732 +- 3.2002

	So we want 0.0573 for Kn and 2.9732 for Ks,
	and the substitution value is between seq 2 and seq 1
	"""
	parts = line.split(" +-")
	ds_value = extract_subs_value(parts[1])
	dn_value = extract_subs_value(parts[0])
	seq_num1 = int(parts[0].split()[0])
	seq_num2 = int(parts[0].split()[1])

	return ds_value, dn_value, seq_num1, seq_num2	

def extract_subs_value(text):
	"""Extract a subsitution value from a line of text.
	
	This is just a friendly function to grab a float value for Ks and Kn 
	values from the junk I get from the last line of the yn00 file.
	
	Line:
	2    1    52.7   193.3   2.0452  0.8979  0.0193 0.0573 +- 0.0177 
	2.9732 +- 3.2002
	
	Parts:
	    ['   2    1    52.7   193.3   2.0452  0.8979  0.0193 0.0573', 
	     ' 0.0177  2.9732', ' 3.2002\n']
	
	So we want 0.0573 for Kn and 2.9732 for Ks.
	"""
	parts = text.split()
	value = float(parts[-1])
	
	return value
	


# nexus works if seprate name and seqs with more than 2 spaces
# phylip doesn't work, has very strange complainments....
# How to prepare a format accepted by yn00??????

class YnCommandline:
	"""Little commandline for yn00.
	"""
	def __init__(self, ctl_file):
		self.ctl_file = ctl_file
		self.parameters = []
		
	def __str__(self):
		 return "yn00 %s" % self.ctl_file
	

