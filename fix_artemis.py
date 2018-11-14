#!/usr/bin/env python
'''
    Fix artemis saved files, with problematic ID missing
/color=([0-9])label=([0-9]+)_(.*)description=(.*)evalue=(.*)hit_id=(.*)$/ID=SS_\2;color=\1;label=\2_\3;description=\4;evalue=\5;hit_id=\6/' 
possible keys:
color, label, evalue, description, alt_score, alt_description, alt_evalue, pseudo_true, part
ID is missing


'''

import sys, re, os
import subprocess
import argparse
from collections import defaultdict
from artemis_string import *


parser = argparse.ArgumentParser()
parser.add_argument('annotation_gff')
parser.add_argument('--prefix', type=str, default="GM_")
args = parser.parse_args()

prefix = args.prefix
annotation_gff = args.annotation_gff
tmp_file = os.path.dirname(os.path.realpath(annotation_gff)) + "/tmp_file"
print tmp_file

def make_gff_attrib(geneid, color, description, pseudo, part):
    label = geneid.split("_")[1] + "_" + description
    attrib = "ID=%s;color=%s;label=%s;description=%s" %(geneid, color, label, description)
    if pseudo:
        attrib += ";pseudo_true=%s" %(pseudo)
    if part:
        attrib += ";part=%s" %(part)
    return attrib


tags = ["color", "label", "evalue", "description", "alt_score", "alt_description", "alt_evalue", "pseudo_true", "part"]

with open(annotation_gff, 'r') as gff_fh, open(tmp_file, 'w') as outh:
    for line in gff_fh:
        line = line.strip()
        arr = line.split("\t")
        source = arr[1]
        feature = arr[2]
        attrib = arr[-1]
        
        if re.match("ID=", attrib):
            print >>outh, line
        else:
            for tag in tags:
                # skip alt_evalue, alt_description if they occur first in the string 
                if tag == "evalue": 
                    i = attrib.find(tag)
                    if attrib[i-4 : i] == "alt_":
                        continue                    
                if tag == "description": 
                    i = attrib.find(tag)
                    if attrib[i-4 : i] == "alt_":
                        # when alt_description is before description
                        attrib = ";description=".join(attrib.rsplit('description=', 1))
                        continue
                    
                attrib = re.sub(tag, ";"+tag, attrib, count=1) # Only the first one, to avoid alt_description
            
            geneid = re.search("label=(\d+)_", attrib).groups()[0]
            geneid = prefix + geneid
            # WB specific
            # if source == "RATT":
            #     geneid = "GL50803_" + geneid
            # elif source == "motif":
            #     geneid = "GI_" + geneid
            # elif source == "PRED":
            #     geneid = "GL_" + geneid
            # else:
            #     print source
                
            attrib = "ID=" + geneid + attrib
            arr[-1] = attrib
            
            print >>outh, "\t".join(arr)

            
    subprocess.call("mv %s %s" %(tmp_file, annotation_gff), shell=True)
