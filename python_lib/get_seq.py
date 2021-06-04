from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC

def get_aa_seq(db_name):
    seq = SeqRecord(Seq(db_name.aa_seq, IUPAC.protein), id=db_name.geneid, description=db_name.description)
    return seq
    
def get_nt_seq(db_name):
    seq = SeqRecord(Seq(db_name.nt_seq, IUPAC.unambiguous_dna), id=db_name.geneid, description=db_name.description)
    return seq


def get_seqs(scf_dict, scfid, start, end, strand, transl_code = 1):
    ''' Retrieve both aa and na seqs from scaffolds dictionary'''
    
    scf_seq = scf_dict[scfid]
    na_seq = scf_seq[start - 1 : end]
    if strand == "-":
        na_seq = na_seq.reverse_complement()
    aa_seq = na_seq.translate(transl_code)
    return na_seq, aa_seq