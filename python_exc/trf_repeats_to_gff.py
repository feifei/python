#!/usr/bin/env python
'''
    Parse TRF html output and convert to gff 
'''
import re, os
import argparse
from collections import defaultdict
from bs4 import BeautifulSoup

def read_soup(infile, n):
    gff_line = ""
    for sub_soup in BeautifulSoup(open(infile), 'lxml'):
        t1 = sub_soup.find('table')
        for r1 in t1.find_all('tr')[1:]:
            pos = r1.find_all('td')[0]
            period = r1.find_all('td')[1]
            copy_num = r1.find_all('td')[2]
            for p, s, c in zip(pos, period, copy_num):
                if re.match("\d", p.get_text()):
                    start, end = map(int, p.get_text().split("--"))
                    n = n + 1
                    attribs = "ID=%s.tr.%d;period_size=%s;copy=%s" %(scfid, n, s.get_text(), c.get_text())
                    gff_line += "%s\t%s\t%s\t%d\t%d\t.\t%s\t.\t%s\n" %(scfid, "trf", "TR", start, end, "+", attribs)
                                    
        break
    return gff_line, n


parser = argparse.ArgumentParser(description='Convert TRF output summary.html file to gff file')
parser.add_argument('summary_html')
parser.add_argument('gff_outfile')
args = parser.parse_args()

summary_html = args.summary_html
gff_outfile = args.gff_outfile
work_dir = os.path.dirname(summary_html)

soup = BeautifulSoup(open(summary_html), "lxml")
table = soup.find('table')
with open(gff_outfile, 'w') as outh:
    for row in table.find_all('tr')[1:]:
        col = row.find_all('td')[1]
        scfid = col.get_text()
        details_file = work_dir + col.a['href']
        n = 0
        gff, n = read_soup(details_file, n)
        print >>outh, gff
        basename, extension = os.path.splitext(details_file)
        new_basename = basename[:-1] + "2"
        extra_file = new_basename + extension
        if os.path.isfile(extra_file):
            gff, n = read_soup(extra_file, n)
            print >>outh, gff

