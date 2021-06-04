'''
    Different overlaps
'''

def overlap(s1, e1, s2, e2):
    overlap_flag = False
    if s1 <= e2 and s2 <= e1:
        overlap_flag = True
    return overlap_flag
    

def complete_overlap(s1, e1, s2, e2):
    overlap_flag = False
    if overlap(s1, e1, s2, e2):
        if (s1 >= s2 and e1 <= e2) or (s1 < s2 and e1 > e2):
            overlap_flag = True
    return overlap_flag


def partial_overlap(s1, e1, s2, e2, p=0.2):
    overlap_flag = False
    if overlap(s1, e1, s2, e2):
        min_size = float(min(e2 - s2, e1 - s1))
        if (s1 >= s2 and s1 <= e2 and (e2 - s1) / min_size >= p) or \
            (s1 <= s2 and s2 <= e1 and (e1 - s2) / min_size >= p):
            # More than p percent mapping
            overlap_flag = True
    return overlap_flag
    