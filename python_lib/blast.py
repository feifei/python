#!/usr/bin/python

import os
import subprocess

def run_formatdb(infile, title, out, args = {}, if_prot='T'):
    """Formatdb. 
    A potential problem is that output folder will be the one runing python """
    # if "data" in os.listdir(os.getcwd()):
    #     os.chdir("data")
    args = ['formatdb', '-i', infile, '-n', out,
            '-p', if_prot, '-t', title] \
            + sum(map(list, zip(args.keys(), args.values())), [])
    cline = subprocess.Popen(args, stdout = subprocess.PIPE)
    cline.stdout.close()
    assert(cline.wait() == 0)
    # os.chdir("..")


def run_blastall(infile, db, program, args = {}, evalue = 0.001, outformat = 0, out= '/dev/stdout'):
    """Blastall """
    args = map(str, ['blastall', '-i', infile, '-p', program, '-d', db,
            '-e', evalue, '-m', outformat, '-o', out] \
            + sum(map(list, zip(args.keys(), args.values())), []))
    cline = subprocess.Popen(args, stdout = subprocess.PIPE)
    cline.stdout.close()
    assert(cline.wait() == 0)


def run_blastp(infile, db, args = {}, evalue = 10, outfmt = 0, out= '/dev/stdout'):
    """Blastall """
    args = map(str, ['blastp', '-query', infile, '-db', db,
            '-evalue', evalue, '-outfmt', outfmt, '-out', out] \
            + sum(map(list, zip(args.keys(), args.values())), []))
    cline = subprocess.Popen(args, stdout = subprocess.PIPE)
    cline.stdout.close()
    assert(cline.wait() == 0)
    