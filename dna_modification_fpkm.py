#!/usr/bin/env python3
'''
    Stats on modification sites
    #m4C/gene, density: # m4C/gene size
    #m6A/gene
    strand vs anti-strand
    #m4C intergenic
    #m6A intergenic
'''

import re, os
import argparse
import gffpandas.gffpandas as gffpd
import pandas as pd
import numpy as np

parser = argparse.ArgumentParser()
parser.add_argument('modification_file')
parser.add_argument('annotation_file')
parser.add_argument('fpkm_file')
parser.add_argument('outfile')

args = parser.parse_args()
annotation_file = args.annotation_file
modification_file = args.modification_file
fpkm_file = args.fpkm_file
outfile = args.outfile
# annotation_file = work_dir + 'app/data/muris.gff'
# modification_file = work_dir + "analysis/modifications/modifications.gff"
# outfile = work_dir + "analysis/modifications/modifications_FPKM.csv"
# fpkm_file = work_dir + 'analysis/rna-seq/rna_cufflinks/genes.fpkm_tracking"


annotation = gffpd.read_gff3(annotation_file).filter_feature_of_type(['gene'])

geneids = (annotation
           .df['attributes']
           .str
           .split(";", n=1, expand=True)[0]
           .str
           .split("=", n=1, expand=True)[1])

description = (annotation
               .df['attributes']
               .str
               .split("product=", n=1, expand=True)[1]
               .str
               .split(";", n=1, expand=True)[0])

pseudo = (annotation
          .df['attributes']
          .str
          .split("pseudo=", n=1, expand=True)[1]
          .str.split(";", n=1, expand=True)[0])

partial = (annotation
           .df['attributes']
           .str
           .split("partial=", n=1, expand=True)[1]
           .str.split(";", n=1, expand=True)[0])


annotation_df = (annotation.df
                 .drop(['attributes'], axis=1)
                 .assign(geneid=geneids)
                 .assign(description=description)
                 .assign(partial=partial)
                 .assign(pseudo=np.where(pseudo.isnull(), 0, 1))) # Make int from pseudo


modifications = gffpd.read_gff3(modification_file).filter_feature_of_type(['m4C', 'm6A'])

def count_modifications(modifications, annotation_df, type='m4C'):
    num_modifications = []
    modifications_subset = modifications.filter_feature_of_type([type])
    for idx, row in annotation_df.iterrows():
        num_modifications.append(
            len(modifications_subset
                .overlaps_with(row['seq_id'], row['start'], row['end'])
                .df))
    return num_modifications

def count_modifications_strand_specific(modifications, annotation_df, type='m4C'):
    num_modifications = []
    modifications_subset = modifications.filter_feature_of_type([type])
    for idx, row in annotation_df.iterrows():
        num_modifications.append(
            len(modifications_subset
                .overlaps_with(row['seq_id'], row['start'], row['end'], strand=row['strand'])
                .df))
    return num_modifications


annotation_modification_df = (annotation_df
                              .assign(m4C=count_modifications(modifications, annotation_df, 'm4C'))
                              .assign(m6A=count_modifications(modifications, annotation_df, 'm6A')))

annotation_modification_df = (annotation_modification_df
                              .assign(m4C_in_strand=count_modifications_strand_specific(modifications, annotation_df, 'm4C'))
                              .assign(m6A_in_strand=count_modifications_strand_specific(modifications, annotation_df, 'm6A')))

annotation_modification_df = (annotation_modification_df
                              .assign(m4C_not_in_strand=annotation_modification_df['m4C'] - annotation_modification_df['m4C_in_strand'])
                              .assign(m6A_not_in_strand=annotation_modification_df['m6A'] - annotation_modification_df['m6A_in_strand']))



# Adding FPKM (from cufflinks) values in the table
fpkm_df = pd.read_table(fpkm_file, sep="\t", header=0)
fpkm_df1 = fpkm_df[['gene_id', 'FPKM']].rename(columns = {'gene_id': 'geneid'})
modifications_fpkm_df = pd.merge(annotation_modification_df, fpkm_df1, on='geneid')

df = (modifications_fpkm_df
      .assign(size=modifications_fpkm_df['end'] - modifications_fpkm_df['start'] + 1))
df = (df
      .assign(m4C_ratio=df['m4C'] / df['size'])
      .assign(m6A_ratio=df['m6A'] / df['size']))


df = df[['geneid', 'pseudo', 'partial', 'seq_id', 'start', 'end', 'size', 'description',
         'm4C', 'm4C_ratio', 'm4C_in_strand', 'm4C_not_in_strand',
         'm6A', 'm6A_ratio', 'm6A_in_strand', 'm6A_not_in_strand',
         'FPKM']]

df.to_csv(outfile)


# pd.set_option('display.max_rows', 500)
# pd.set_option('display.max_columns', 100)
# pd.set_option('display.width', 1000)


print("m4C : %d" %len(modifications.df.loc[modifications.df['type'] == 'm4C']))
print("m6A : %d" %len(modifications.df.loc[modifications.df['type'] == 'm6A']))









