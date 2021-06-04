#!/usr/bin/python

from Bio.Seq import Seq

def translate_longest(sequence, table="Standard", stop_symbol="*", to_stop=False):
    sequence = sequence[0: len(sequence) - len(sequence) % 3]
    aa = sequence.translate(table, stop_symbol, to_stop)
    return sorted(aa.split('*'), lambda x,y: cmp(len(y), len(x)))[0]


def translate_longest_from_M(sequence, table="Standard", stop_symbol="*", to_stop=False):
    sequence = sequence[0: len(sequence) - len(sequence) % 3]
    aa = sequence.translate(table, stop_symbol, to_stop)
    return sorted(aa.split('*'), lambda x,y: cmp(len(y[y.find("M"):]), len(x[x.find("M"):])))[0]
    # max_aa = sorted(aa.split('*'), lambda x,y: cmp(len(y[y.find("M"):]), len(x[x.find("M"):])))[0]
    # return max_aa[max_aa.find("M"): ]


def longest_orf(sequence, table="Standard", stop_symbol="*", from_M = False, to_stop=False):
    sequence = sequence[0: len(sequence) - len(sequence) % 3]
    aa = sequence.translate(table, stop_symbol, to_stop)
    if from_M:
        aa_l = translate_longest_from_M(sequence, table, stop_symbol, to_stop)
    else:
        aa_l = translate_longest(sequence, table, stop_symbol, to_stop)
    start = aa.find(aa_l)
    end = start + len(aa_l) + 1 # + 1 to include stop codon
    start, end = 3 * start + 1, 3 * end
    return (sequence[start - 1: end], start, end) 




def longest_orf_nearby(sequence, ref_length, qry_start, qry_end, table="Standard", stop_symbol="*", start_symbol="M", to_stop=False):
    ''' Finding the longest orf extending the given positions on the given frame
        ref_length in aa
    '''
    aa = sequence.translate(table, stop_symbol, to_stop)
    start, end = qry_start / 3 + 1, qry_end / 3
    if aa[start - 1 : end - 1].find(stop_symbol) > 0:
        return None
    else:
        if aa[end : ].find(stop_symbol) < 0:
            end = len(aa)
        for i, c in enumerate(aa[end - 1 : ]):
            if c == stop_symbol:
                end = end + i
                break
        
        starts = [] # Include the first codon in the given stretch
        for i, c in enumerate(aa[0 : start][::-1]):
            if c == stop_symbol:
                starts.append(start - i)
                break
            elif c == start_symbol:
                s = start - i - 1
                starts.append(s)
                if (end - s) > ref_length * 2:
                    break
                    
        # Pick the start which gives the cloest length to ref_length
        best_length = 0
        for s in starts:
            aa_length = end - start + 1
            if abs(aa_length - ref_length) < abs(best_length - ref_length):
                best_length = aa_length
                start = s
        
        aa_seq = aa[start: end]
        if aa_seq and aa_seq.count("X") / float(len(aa_seq)) < 0.2 and len(aa_seq) > 50: 
            #aa_seq.endswith(stop_symbol) and aa_seq.startswith(start_symbol) and abs(len(aa_seq) - ref_length) / float(ref_length) < 0.2:
            start, end = start * 3 + 1, end * 3
            na_seq = sequence[start - 1 : end]
            return (na_seq, aa_seq, start, end)
        else:
            return None

    
Seq.translate_l = translate_longest
Seq.translate_l_from_M = translate_longest_from_M
Seq.longest_orf = longest_orf
Seq.longest_orf_nearby = longest_orf_nearby

