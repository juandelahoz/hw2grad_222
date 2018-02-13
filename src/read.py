import numpy as np

class read():

    read_id  = None
    read_len = None
    sequence = None

    aligned  = False
    position = "."
    strand   = "."
    mismatch = ["."]
    deletion = ["."]
    insert   = ["."]


    pair_position = "."
    pair_strand   = "."

    def __init__(self, rdid, seq):
        self.sequence = seq
        self.read_id = rdid
        self.read_len = len(seq)

# TODO: update summary
    def summary(self):
        """returns a string with a summary of the alignment of the read"""
    # for deletions
        dels_p = ["."]
        dels_a = ["."]
        if deletion != ["."]:
            deletions = deletion.keys()
            for dtn in deletions:
                dels_p.append(      str(dtn) )
                dels_a.append( deletion[dtn] )
    # for insertions
        ins_p = ["."]
        ins_a = ["."]
        if insert != ["."]:
            insertion = insert.keys()
            for itn in insertion:
                ins_p.append(      str(itn) )
                ins_a.append( deletion[itn] )


        return '{}\t{}\t{}\t{}\t{}\t{}\t{}'.format(self.sequence,self.read_id,
            self.strand,self.position,
            ",".join(map(str, self.mismatch)),
            ",".join(map(str, self.dels_p)),",".join(map(str, self.dels_a)),
            ",".join(map(str, self.ins_p )),",".join(map(str, self.ins_a )),
            self.pair_strand,self.pair_position)

