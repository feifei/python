#!/usr/bin/python

import re
from itertools import tee, islice, chain, izip
from difflib import SequenceMatcher


# convert string to bool
def str2bool(s):
    return bool(eval(s.capitalize()))


# Strip the last extension
def strip_extension(str):
   return re.sub('\.[^.]+$', '', str)


# Strip the file names, leaving the directory
def strip_filename(str):
   return re.sub('/[^/]+$', '', str)


# return a list of region where the values lie inside the give edge values
def find_region(L, edge):
	if len(edge) != 2:
		small = big = edge[0]
	else:
		small = edge[0]
		big = edge[1]
	start = 0
	stop = 0
	if big >= L[-1]:
		stop = len(L)-1
	for i, l in enumerate(L):
		if l <= small:
			start = i
		elif l > big:
			stop = i
			break

	return range(start, stop+1)
	#return range(start, stop+1)

# return a uniq list of the give list
def uniq_list(seq):  
	# order preserving 
	checked = [] 
	for e in seq: 
		if e not in checked: 
		 	checked.append(e) 
	return checked

# order given list by its give index list
def order_index(L, index):
	new_L = []
	for i in index:
		new_L.append(L[i])

	return new_L
    

def atoi(text):
    return int(text) if text.isdigit() else text

def natural_keys(text):
    '''
    alist.sort(key=natural_keys) sorts in human order
    http://nedbatchelder.com/blog/200712/human_sorting.html
    (See Toothy's implementation in the comments)
    '''
    return [ atoi(c) for c in re.split('(\d+)', text) ]


def previous_and_next(some_iterable):
    prevs, items, nexts = tee(some_iterable, 3)
    prevs = chain([None], prevs)
    nexts = chain(islice(nexts, 1, None), [None])
    return izip(prevs, items, nexts)



def similar(a, b):
    return SequenceMatcher(None, a, b).ratio()


def gene_to_scf_pos(strand, scf_begin, scf_end, gene_begin, gene_end):
    # Convert nt position at gene level to position at scaffold level
    if strand == "+":
        begin_global = scf_begin + gene_begin - 1
        end_global = scf_begin + gene_end - 1
    elif strand == "-":
        begin_global = scf_end - gene_end + 1
        end_global = scf_end - gene_begin + 1
    return begin_global, end_global


def aa_to_na_pos(begin, end):
    # Convert position at amino acid level to nucleotide level
    # 1-based
    return begin * 3 - 2, end * 3 