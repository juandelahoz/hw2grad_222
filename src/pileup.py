"""Script to generate a pileup from aligned reads and a reference"""
import sys
import util    as u

def pileup(aligns, ref_genome):
    """take a file with aligned reads and extract variants."""
    print("Generating Pileup from file: ")
    print("   " + aligns)

    # open aligned reads file, and pileup file
    reads_f  = open(aligns, 'r')
    pileup_f = open(aligns + ".plp", "w")

    # set variables
    genome = ref_genome.getSequence()[0]
    pile   = [0] * len(genome)
    sym    = {"A":0,"C":1,"G":2,"T":3}
    let    = {0:"A",1:"C",2:"G",3:"T"}
    donor  = [list(pile),list(pile),list(pile),list(pile)]
    end    = False         # flag for end of reads
    window = 1000          # size of the window to slide

    # write with the expected format
    pileup_f.write(">" + list(ref_genome.genome.keys())[0] +"\n")
    pileup_f.write(">SNP\n")

    start = 0
    reads = []
    # keep function running while there are reads in file
    while not end:

        # collect the information of all the reads inside the window
        for read in reads_f:
            read = read.strip().split('\t')
            seq,init,pos = read[0],int(read[3]),read[2]
            if pos == "-":
                seq = u.rev(seq)                   # reverse if inverted
            reads.append((seq,init))               # save the read
            if init >= start + window:             # stop with the first that scapes
                break                              #           the end of the window

        for read in reads:
            seq  = read[0]                         # retrieve the reads from array
            init = read[1]
            for i in range(len(seq)):              # position by position of each read
                donor[ sym[seq[i]] ][init+i] += 1  #  add to PILEUP!
        
        # call variants:
        depth = 0
        dom = 0
        final = min(len(genome), start + window)
        for j in range( start, final ):     # for each position in the window ...
            don = genome[j]                      # by default, the donor is = to reference
            for nt in range(4):                  # for each nucleotide... (ACGT)
                depth += donor[nt][j]        # add to sequencing depth for position j
                if dom < donor[nt][j]:
                    dom = donor[nt][j]   # save the max depth for any nucleotide in pos j
                    don = let[nt]        #  and the identity of that nucleotide
            if don != genome[j]:         # if donor and reference diverge...
                if depth > 2:                                # QC: and donor has >2 reads
                    if float(dom)/float(depth) > 0.6:        # QC: and freq(nt) >0.6
                        #print(genome[j],donor[0][j],donor[1][j],donor[2][j],donor[3][j],don,j)
                        pileup_f.write(genome[j] +","+ don +","+ str(j) +"\n")
            depth = 0
            dom = 0
        
        # reset variables or finish outputing
        if len(genome) > start + window:
            start += window
            reads = []
        else:
            end = True

    pileup_f.close()
    return 1
