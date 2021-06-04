from Bio.Seq import Seq
from translate_l import *

def normal_pos_to_frame(seq, f, p):
    if f < 0:
        return len(seq) + f - p + 2
    else:
        return p - (f - 1)

def frame_pos_to_normal(seq, f, p):
    if f < 0:
        return len(seq) + f - p + 2
    else:
        return p + (f - 1)

def frame(seq, f):
    if f < 0:
        return seq.reverse_complement()[-1 - f:]
    else:
        return seq[f - 1:]

def frames(seq, fs = [1,2,3,-1,-2,-3]):
    return map(lambda f: frame(seq, f), fs)

def best_frame(seq, trans_code = 6):
    best_f = 0
    best_len = 0
    for f in [1,2,3,-1,-2,-3]:
        nt = seq.frame(f)
        aa_l = nt.translate_l(trans_code)
        if len(aa_l) > best_len:
            best_len = len(aa_l)
            best_f = f
    return best_f


def best_frame_from_M(seq, trans_code = 6):
    best_f = 0
    best_len = 0
    for f in [1,2,3,-1,-2,-3]:
        nt = seq.frame(f)
        aa_l = nt.translate_l_from_M(trans_code)
        if len(aa_l) > best_len:
            best_len = len(aa_l)
            best_f = f
    return best_f


Seq.normal_pos_to_frame = normal_pos_to_frame
Seq.frame_pos_to_normal = frame_pos_to_normal
Seq.frame = frame
Seq.frames = frames
Seq.best_frame = best_frame
Seq.best_frame_from_M = best_frame_from_M

