import pysam

# Return if a read is clipped on the left
def is_left_clipped(read):
    if read.cigartuples is not None:
        if read.cigartuples[0][0] in [4, 5]:
            return True
    return False


# Return if a read is clipped on the right
def is_right_clipped(read):
    if read.cigartuples is not None:
        if read.cigartuples[-1][0] in [4, 5]:
            return True
    return False


# Return if a read is clipped on the right or on the left
def is_clipped(read):
    if read.cigartuples is not None:
        if is_left_clipped(read) or is_right_clipped(read):
            return True
    return False


# Return the mate of a read. Get read mate from BAM file
def get_read_mate(read, bamfile):
    # print(read)
    # print('Mate is mapped')
    iter = bamfile.fetch(read.next_reference_name,
                         read.next_reference_start, read.next_reference_start + 1,
                         multiple_iterators=True)
    for mate in iter:
        if mate.query_name == read.query_name:
            if (read.is_read1 and mate.is_read2) or (read.is_read2 and mate.is_read1):
                # print('Mate is: ' + str(mate))
                return mate
    return None