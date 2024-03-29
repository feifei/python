#!/usr/bin/env python
'''
assemstats.py

Uses screed to calculate various assembly statistics.

You can obtain screed at github by running
   git clone git://github.com/acr/screed.git

Then, install by running
   python setup.py install
in the newly created screed directory.

Once completed, you should be able to run this script as is.

Author: Jason Pell (pelljaso@cse.msu.edu)
'''

import screed.fasta
import sys, re
import glob

def trimLens(lens, minLen):
   '''
   Eliminates any reads below a certain threshold.  Function assumes that input
   list lens is sorted smallest to largest.
   '''

   index = 0

   for i in range(len(lens)):
      if lens[i] < minLen:
         index += 1
      else:
         break

   return lens[index:len(lens)]

def getGC(seq):
    return sum([1.0 for nucl in seq if nucl in ['G', 'C', 'g', 'c']])

def getLens(filename):
   '''
   Parses FASTA file using screed to create a sorted list of contig lengths.
   '''
   lens = []
   lens_without_N = []
   gcs = []
   fd = open(filename, 'r')

   fa_instance = screed.fasta.fasta_iter(fd)
   for record in fa_instance:
      lens.append(len(record['sequence']))
      lens_without_N.append(len(re.sub("N+", "", record['sequence'].upper())))
      gcs.append(getGC(record['sequence']))

   fd.close()

   return sorted(lens), sorted(lens_without_N), gcs

def calcNXX(lens, percent):
   '''
   Calculates any NXX (e.g. N50, N90) statistic.
   '''

   lenSum = sum(lens)
   threshold = (float(percent) / 100) * lenSum
   
   runningSum = 0
   nxx = 0
   nxxLen = 0

   for i in range(len(lens)-1, -1, -1):
      myLen = lens[i]
      nxx += 1
      runningSum += myLen

      if runningSum >= threshold:
         nxxLen = myLen
         break

   return nxx, nxxLen



def main():
   '''
   Outputs assembly statistics for provided FASTA files.
   '''
   if len(sys.argv) < 3:
      print "Usage: python assemstats.py <min contig length> [ FASTA files ]"
      return

   try:
      minLen = int(sys.argv[1])
   except ValueError:
      print "Minimum contig length must be an integer."
      return

   #print "filename sum n trim_n min med mean max n50 n50_len n90 n90_len"

   for expr in sys.argv[2:len(sys.argv)]:
      for filename in glob.glob(expr):
         lens, lens_without_N, gcs = getLens(filename)
         trimmedLens = trimLens(lens, minLen)
         trimmedLens_without_N = trimLens(lens_without_N, minLen)

         if len(trimmedLens) == 0:
            print filename + " - no sequences longer than threshold"
            continue

         statN = len(lens)
         statTrimmedN = len(trimmedLens)
         statSum = sum(trimmedLens)
         statSum_without_N = sum(trimmedLens_without_N)
         statMin = min(trimmedLens)
         statMax = max(trimmedLens)
         statMed = trimmedLens[(len(trimmedLens)-1)/2]
         statMean = int(statSum / float(statTrimmedN))
         statN50, statN50Len = calcNXX(trimmedLens, 50)
         statN75, statN75Len = calcNXX(trimmedLens, 75)
         statN80, statN80Len = calcNXX(trimmedLens, 80)
         statN90, statN90Len = calcNXX(trimmedLens, 90)
 
         print str(filename)
         print "Total number of contigs: %d" %statN
         print "Total number of contigs >= %d: %d" %(minLen, statTrimmedN)
         print "Total bp after trimming: %d" %statSum
         print "Total bp after trimming without N: %d" %statSum_without_N
         print "Minumum contig length: %d" %statMin
         print "Median contig length: %d" %statMed
         print "Mean contig length: %d" %statMean
         print "Biggest contig length: %d" %statMax
         print "N50 inclues %d contigs and N50 length: %d" %(statN50, statN50Len)
         print "N75 inclues %d contigs and N75 length: %d" %(statN75, statN75Len)
         print "N80 inclues %d contigs and N80 length: %d" %(statN80, statN80Len)
         print "N90 includes %d contigs and N90 length: %d" %(statN90, statN90Len)
         print "GC: %.1f" %(sum(gcs)/sum(lens) * 100) 

main()