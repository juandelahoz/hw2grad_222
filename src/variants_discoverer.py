"""Script to find variants from a donor in NGS paired end reads."""
import sys
import numpy   as np
import util    as u
import read    as rd
import genome  as gm
import aligner as aln
import pileup  as plp

def main():
    """Main entry point for the script."""
    ref_fa = sys.argv[1]
    reads  = sys.argv[2]
    ref_seq = load_reference(ref_fa)

#    align_reads(reads, ref_seq) 
    plp.pileup(reads, ref_seq)
#    ref_seq.index_BWT(2000)
    

# TODO: add option for unpaired reads
def align_reads(read_fn, ref_genome):
    """take a file with reads and align read by read."""
    print("Aligning reads from file: ")
    print("   " + read_fn)

    # open "fastq" file to align, and "bam" to output
    fastq   = open(read_fn, 'r')
    mapped  = open(read_fn + ".aln", "w")

    read_id = 0
    for read in fastq:
        if read_id > 0:              # Skip the first line

            # read, separate and prepare reads
            read = read.strip()
            reads_raw = read.split(',')
            read_1 = rd.read(str(read_id)+"_1", reads_raw[0])
            read_2 = rd.read(str(read_id)+"_2", reads_raw[1])

            # align and output
            aln_reads = aln.align_trivial([read_1,read_2], ref_genome)
            for alignment in aln_reads:
                mapped.write( alignment.summary() +"\n")

        if read_id % 10 == 0:              # report progress
            print("{} reads aligned".format(read_id))
        read_id += 1

    mapped.close()
    return 1

def load_reference(ref_fn):
    fasta = open(ref_fn, 'r')

    # set variables:
    seqfa = []
    seqid = []
    seq = ""
    fline = False

    for line in fasta:
        if line.startswith(">"):              # for lines with ids...
            seqid.append(line.strip()[1:])    #   save id
            if fline:
                seqfa.append(seq)             #   save seq and
                seq = ""                      #     start a new seq
            fline = True
        else:                                 # for lines with sequence...
            seq += line.strip()
    seqfa.append(seq)

    # initialize the genome object
    genome = gm.genome(seqid[0],seqfa[0])

    # add other sequences, if any
    for i in range(1,len(seqid)):
        genome = gm.genome.add_sequence(seqid[i],seqfa[i])

    # report progress
    print("Reference genome loaded.")
    print("         total length: {}".format(genome.length))
    return genome

def print_variants():
    """format and print the list of variatns"""
    pass

if __name__ == '__main__':
    sys.exit(main())
