"""Script to generate a pileup from aligned reads and a reference"""
import sys
import util    as u

def pileup(aligns, ref_genome):
    """take a file with aligned reads and extract variants."""
    print("Generating Pileup from file: ")
    print("   " + aligns)

    # open aligned reads file, and pileup file
    reads_f  = open(aligns, 'r')

    # set variables
    genome = ref_genome.getSequence()[0]
    pile   = [0] * len(genome)

    start = 0
    reads = []
    window = 1000          # size of the window to slide

    snps = {}
    dels = {}              # format: key=position, vals=(ref, alt, support)
    inse = {}

    # collect the information of all the reads inside the window
    for read in reads_f:
        read = read.strip().split('\t')
        seq,init,pos = read[0],int(read[3]),read[2]
        if pos == "-":
            seq = u.rev(seq)                   # reverse if inverted

        # get differences if any:
        # SNPs
        if read[4] != ".":
            for i in map(int, read[4].split(",")):
                if (init + i) in snps:        # add support if it exists
                        # pos_on_ref      ref_allele  alt_allele         support +1
                    snps[init + i] = (genome[init + i],   seq[i],   (snps[init+i][2]+1) )
                else:                          # or create the event if it doesn't
                    snps[init + i] = (genome[init + i],   seq[i], 1)
        # deletions
        if read[5] != ".":
            pass
        # insertions
        if read[7] != ".":
            pass

        # write with the expected format
    pileup_f = open(aligns + ".plp", "w")

    pileup_f.write(">" + list(ref_genome.genome.keys())[0] +"\n")
    pileup_f.write(">SNP\n")
    snps_pos = list(snps.keys())
    snps_pos.sort()
    for j in snps_pos:
        if snps[j][2] > 3:
            pileup_f.write( snps[j][0] + "," + snps[j][1] + "," + str(j) + "\n")

    pileup_f.close()

    return 1
