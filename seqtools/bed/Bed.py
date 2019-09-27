def count_bed(bed, strand=None):
    '''Counts number of entry in BED, can be limited to a specific strand.'''
    count = 0
    with open(bed, "r") as infile:
        for line in infile:
            if line.startswith('track') or line.startswith('browser') or line.startswith('#'):
                continue
            if strand is None:
                count += 1
            else:
                columns = line.rstrip('\r\n').split('\t')
                if len(columns) >= 6 and columns[5] == strand:
                    count += 1
    return count


def empty_bed(bed_output, sample, strand=None):
    '''Create an empty BED file.'''
    track = 'track type=bedGraph name="' + sample + '"'
    if not strand is None:
        track += ' Minus' if strand == '-' else ' Plus'
    with open(bed_output, "w") as outfile:
        outfile.write(track + '\n')
