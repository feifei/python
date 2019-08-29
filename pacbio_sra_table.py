#!/usr/bin/env python
'''
    Create table for the submission of SRA data
'''

import os, re
import argparse
import glob

parser = argparse.ArgumentParser(description='Process .h5 files to SRA table')
parser.add_argument('--bioproject', dest = "bioproject", default = "PRJNA531391")
parser.add_argument('--biosample', dest = "biosample", default = "SAMN02952905")
parser.add_argument('--title', dest = "title", default = "WB_DNA")
parser.add_argument('--design_description', dest = "design_description", default = "giardia intestinalis WB DNA WGS")
parser.add_argument('--dir', dest='dir', default = "/Volumes/Molev/data/pacbio/pb_381/rawdata",
                    help='specify dir where h5 files are')
args = parser.parse_args()
bioproject = args.bioproject
biosample = args.biosample
title = args.title
design_description = args.design_description
dir = args.dir
#h5_files = "/Volumes/Molev/data/pacbio/pb_381/rawdata/*/*/*.h5"
#listing = glob.glob(h5_files)
study_id = re.match(".*?/pacbio/(.*?)/rawdata", dir).groups()[0]
for _, dirs, _ in os.walk(dir):
    dirs = sorted(dirs)
    for subdir in dirs:
        cell_id = subdir
        library_id = study_id + "_" + cell_id
        h5_dir = dir + "/" + subdir + "/Analysis_Results/"
        filenames = ""
        for _, _, files in os.walk(h5_dir):
            filenames = "\t".join(sorted(files))                
        print "\t".join([bioproject, biosample, library_id, title, "WGS", "GENOMIC", "PCR", \
        "single", "PACBIO_SMRT", "PacBio RS II", design_description, "PacBio_HDF5", filenames])
    
    break
