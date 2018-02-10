
def test_read_format(read):
    """get reads file and check if their format is ok"""
    pass

def revCompl(read):
    comp = {"A":"T", "T":"A", "C":"G", "G":"C"}
    rev_comp = ""
    for nucl in read:
        rev_comp = comp[nucl] + rev_comp
    return rev_comp

def rev(read):
    rev = ""
    for nucl in read:
        rev = nucl + rev
    return rev

def SymbolToNumber(symbol):
    symbols = {"A":0,"C":1,"G":2,"T":3}
    return symbols[symbol]

def NumberToSymbol(number):
    numbers = {0:"A",1:"C",2:"G",3:"T"}
    return numbers[number]

def PatternToNumber(pattern):
    if pattern == "":
        return 0
    symbol = pattern[-1]
    prefix = pattern[:-1]
    return 4 * PatternToNumber(prefix) + SymbolToNumber(symbol)

def NumberToPattern(index, k):
    if k == 1:
        return NumberToSymbol(index)
    prefixIndex = int(index/4)
    r = index%4
    symbol = NumberToSymbol(r)
    PrefixPattern = NumberToPattern(prefixIndex,k-1)
    return PrefixPattern + symbol



def unpermute_BWT(index,count):
    # try to reconstruct the original sequence (just to test if the index is correct!)
    # read the index
    nt = ["A","C","G","T"]
    # start unpermuting from the first position
    seq = index[0][0]
    pos = index[0][1]
    while seq[0] != "$":         # until the first character
        cpos = 0
        for n in nt:
            if n != seq[0]:
                cpos += count[n] # advance first column index per nt
            else:
                break
        cpos += pos              # advance to the position of the nt
        seq = index[cpos][0] + seq
        if index[cpos][0] != "$":
            pos = index[cpos][1]
    return seq[1:]               # omit initial '$'

def pretty_print_aligned_reads_with_ref(genome_oriented_reads, read_alignments, ref, read_length=50,
                                        line_length=100, read_sep=100, buffer=30):
    """
    :param genome_oriented_reads: oriented reads generated by an alignment algorithm
    :param read_alignments: alignments generated from an alignment algorithm
    :param ref: reference generated by read_ref
    :return: Returns nothing, but prints the reads aligned to the genome to
     show you what pileup actually *LOOKS* like. You should be able to call SNPs
     by eyeballing the output. However, there are some reads that will not align.
     In the future you'll want to re-check why these reads aren't aligning--the cause
     is usually a structural variation, like an insertion or deletion.
    """
    output_str = ''
    good_alignments = [read_sep + read_length - buffer < x[1] - x[0] <
                       read_sep + read_length + buffer for x in read_alignments]
    # There should be read_length + x (90 < x < 110) p between the reads, and we give a little
    # extra space in case there's been a deletion or insertion.  Depending on the type of
    # deletions/insertions

    best_reads = [genome_oriented_reads[i] for i in range(len(good_alignments))
                  if good_alignments[i]]
    # Remove the reads that do not have a good alignment, or a good reverse alignment.
    best_alignments = [read_alignments[i] for i in range(len(read_alignments))
                       if good_alignments[i]]
    # Take their corresponding alignments
    aligned_reads = [best_reads[i][0] + '.' * (best_alignments[i][1] - best_alignments[i][0] - read_length)
                     + best_reads[i][1] for i in range(len(best_reads))]
    # This turns the reads into strings oriented towards the genome.
    # We get the first read, followed by the correct number of dots to join the first and second reads,
    # and then the second read.

    first_alignment = [x[0] for x in best_alignments]
    alignment_indices = np.argsort(first_alignment)
    sorted_reads = np.array([aligned_reads[i] for i in alignment_indices])
    sorted_alignments = np.array([best_alignments[i] for i in alignment_indices])

    # You don't need to worry too much about how the code block below works--its job is to make it so
    # that a read that starts printing in the third row will continue printing in the third row of the
    # next set of lines.
    active_reads = []
    output_str += '\n\n' + '-' * (line_length + 6) + '\n\n'
    read_indices = np.array([sorted_alignments[j][0]/line_length for j in range(len(sorted_alignments))])

    for i in range(len(ref) / line_length):
        next_ref = ref[i * line_length: (i + 1) * line_length]
        read_mask = (read_indices == i)
        new_alignments = sorted_alignments[read_mask]
        new_reads = sorted_reads[read_mask]
        space_amounts = [_[0] % line_length for _ in new_alignments]
        new_reads_with_spaces = [' ' * space_amounts[j] + new_reads[j] for j in range(len(new_reads))]
        empty_active_read_indices = [index for index in range(len(active_reads)) if active_reads[index] == '']
        for j in range(min(len(new_reads_with_spaces), len(empty_active_read_indices))):
            active_reads[empty_active_read_indices[j]] = new_reads_with_spaces[j]

        if len(new_reads_with_spaces) > len(empty_active_read_indices):
            active_reads += new_reads_with_spaces[len(empty_active_read_indices):]
        printed_reads = ['Read: ' + read[:line_length] for read in active_reads]
        active_reads = [read[line_length:] for read in active_reads]
        while len(active_reads) > 0:
            last_thing = active_reads.pop()
            if last_thing != '':
                active_reads.append(last_thing)
                break
        output_lines = ['Ref:  ' + next_ref] + printed_reads
        output_str += 'Reference index: ' + str(i * line_length) + \
                      '\n' + '\n'.join(output_lines) + '\n\n' + '-' * (line_length + 6) + '\n\n'
    # print output_str
    return output_str
