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
    snps = {}
    dels = {}   # format: key=position, vals=(ref, (alt,) support)
    insr = {}

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
                # add support if it exists (position and allele)
                if (init + i) in snps:
                    if seq[i] == snps[init + i][1]:
                            # pos_on_ref      ref_allele  alt_allele         support +1
                        snps[init + i] = (genome[init + i],   seq[i],   (snps[init+i][2]+1) )
                else:                          # or create the event if it doesn't
                    snps[init + i] = (genome[init + i],   seq[i], 1)

        # deletions
        if read[5] != ".":
            pos_d = list(map(int, read[5].split(",")))
            ale_d = list(map(int, read[6].split(",")))
            for i in range(len(pos_d)):
                pi_del = init + pos_d[i] + 1
                pf_del = init + pos_d[i] + 1 + ale_d[i]
                # add support if it exists (position and allele)
                if pi_del in dels:
                    if dels[pi_del][0] == genome[ pi_del:pf_del ]:
                        dels [pi_del][1] += 1
                else:                          # or create the event if it doesn't
                    dels[pi_del] = [genome[ pi_del:pf_del ],  1]

        # insertions
        if read[7] != ".":
            pos_i = list(map(int, read[7].split(",")))
            ale_i = list(map(int, read[8].split(",")))
            for i in range(len(pos_i)):
                pi_ins = pos_i[i] + 1
                pf_ins = pos_i[i] + 1 + ale_i[i]
                p_geno = pos_i[i] + 1 + init
                # add support if it exists (position and allele)
                if p_geno in insr:
                    if insr[p_geno][0] == seq[ pi_ins:pf_ins ]:
                        insr[p_geno][1] += 1
                else:                          # or create the event if it doesn't
                    insr[p_geno] = [seq[ pi_ins:pf_ins ],  1]

        # write with the expected format
    pileup_f = open(aligns + ".plp", "w")
    pileup_f.write(">" + list(ref_genome.genome.keys())[0] +"\n")

    pileup_f.write(">INS\n")
    insr_pos = list(insr.keys())
    insr_pos.sort()
    for j in insr_pos:
        if insr[j][1] > 2:
            pileup_f.write( insr[j][0] + "," + str(j) + "\n")

    pileup_f.write(">DEL\n")
    dels_pos = list(dels.keys())
    dels_pos.sort()
    for j in dels_pos:
        if dels[j][1] > 2:
            pileup_f.write( dels[j][0] + "," + str(j) + "\n")

    pileup_f.write(">SNP\n")
    snps_pos = list(snps.keys())
    snps_pos.sort()
    for j in snps_pos:
        if snps[j][2] > 3:
            pileup_f.write( snps[j][0] + "," + snps[j][1] + "," + str(j) + "\n")

    pileup_f.close()

    return 1
