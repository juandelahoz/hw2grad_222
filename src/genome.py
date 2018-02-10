import util as u

class genome():

    length = 0
    kmers = None
    genome = {}


    def __init__(self, ids, seqs):
        self.genome[ids] = seqs
        self.length = len(seqs)
        self.index_kmers()

    def add_sequence(ids, seqs):
        self.genome[ids] = seqs
        self.length += len(seqs)

    def index_kmers(self):
        """creates the index of the reference"""
        ref = self.genome[list(self.genome.keys())[0]]
        k = 6        # index positionos of kmers
        index = [None] * (4**k)
        for i in range(len(ref)-k+1):
            kmer_num = u.PatternToNumber( ref[i:i+k] )
            if(index[kmer_num] == None):
                index[kmer_num] = [i]
            else:
                index[kmer_num].append(i)
        self.kmers = index

    def getSequence(self):
    	chrs = list(self.genome.keys())
    	sequence = []
    	for chrom in chrs:
    		sequence.append(self.genome[chrom])
    	return sequence

    def index_BWT(self,k): # include the kmer size
        fastaidx = open(list(self.genome.keys())[0] + ".idx", "w")

        sequence = self.getSequence()
        sequence = sequence[0] + "$"  # for now its just a single chrom
        L = len(sequence) # save the length of reference

        # generate the index from kmers
        index_m = []
        for i in range(L):
            if i < L-k+1:
                index_m.append((sequence[i:i+k],i))
            else:
                index_m.append((sequence[i:] + sequence[:i-(L-k)],i))

        index_m.sort()      # this is the critical step!

        index = ""
        for kmer in index_m:
            if kmer[1] == 0:
                index += "$"   # remember that before the first is the last
            else:
                index += sequence[kmer[1]-1]

        s = 1               # keep track of the position of letters each "s"
        count = {"A":0,"C":0,"G":0,"T":0}
        for nt in index:
            if nt != "$":
                count[nt] += 1
                if count[nt] % s == 0:    # report counts every "s" nt only
                    fastaidx.write('{} {}\n'.format(nt,str(count[nt])))
                else:                
                    fastaidx.write(nt + "\n")
            else:
                fastaidx.write(nt + "\n")
