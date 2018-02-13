class read():

    read_id  = None
    read_len = None
    sequence = None

    aligned   = False
    position  =  "."
    strand    =  "."
    mismatch  = ["."]
    deletion  = ["."]
    insertion = ["."]

    pair_position = "."
    pair_strand   = "."

    def __init__(self, rdid, seq):
        self.sequence = seq
        self.read_id = rdid
        self.read_len = len(seq)

    def summary(self):
        """returns a string with a summary of the alignment of the read"""
    # for deletions
        dels_p = ["."]
        dels_a = ["."]
        if self.deletion != ["."]:
            deletions = self.deletion.keys()
            for dtn in deletions:
                dels_p.append(           str(dtn) )
                dels_a.append( self.deletion[dtn] )
    # for insertions
        ins_p = ["."]
        ins_a = ["."]
        if self.insertion != ["."]:
            insertions = self.insertion.keys()
            for itn in insertions:
                ins_p.append(            str(itn) )
                ins_a.append( self.insertion[itn] )

        return '{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}'.format(
            self.sequence,self.read_id,
            self.strand,self.position,
            ",".join(map(str, self.mismatch)),
            ",".join(map(str, dels_p)),",".join(map(str, dels_a)),
            ",".join(map(str, ins_p )),",".join(map(str, ins_a )),
            self.pair_strand,self.pair_position)
