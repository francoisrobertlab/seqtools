import logging
import os
import subprocess

import click

BASE_SCALE = 1000000


@click.command()
@click.option('--samples', '-s', type=click.File('r'), default='samples.txt',
              help='Sample names listed one sample name by line.'
              ' The first line is ignored.'
              ' An SRR id can be provided (tab-separated) to download the FASTQ file automatically, otherwise it must be provided.')
@click.option('--fasta', '-f', type=click.Path(exists=True), default='sacCer3.fa',
              help='FASTA file used for alignment.')
@click.option('--sizes', type=click.Path(exists=True), default='sacCer3.chrom.sizes',
              help='Size of chromosomes.')
@click.option('--threads', '-t', default=1,
              help='Number of threads used to process data.')
def main(samples, fasta, sizes, threads):
    '''Analyse Martin et al. data from November 2018 in Genetics.'''
    logging.basicConfig(filename='debug.log', level=logging.DEBUG)
    bwa_index(fasta)
    next(samples)
    samples_lines = samples.read().splitlines()
    for sample_line in samples_lines:
        if sample_line.startswith('#'):
            continue
        sample_info = sample_line.split('\t');
        sample = sample_info[0]
        srr = sample_info[1] if len(sample_info) > 1 else None
        analyse(sample, srr, fasta, sizes, threads)


def analyse(sample, srr, fasta, sizes, threads):
    '''Analyse a single sample.'''
    print ('Analyse sample {}'.format(sample))
    fastq = sample + '_1.fastq'
    fastq2 = sample + '_2.fastq'
    if not os.path.isfile(fastq):
        fastq = sample + '_1.fastq.gz'
        fastq2 = sample + '_2.fastq.gz'
    if not os.path.isfile(fastq):
        fastq = sample + '.fastq'
    if not os.path.isfile(fastq):
        fastq = sample + '.fastq.gz'
    if not os.path.isfile(fastq) and srr:
        fastq = sample + '_1.fastq.gz'
        fastq2 = sample + '_2.fastq.gz'
        print ('Downloading FASTQ for {} with SRR {}'.format(sample, srr))
        download(sample, srr)
    bam_raw = sample + '-raw.bam'
    print ('Running BWA on sample {}'.format(sample))
    bwa(fastq, fastq2, fasta, bam_raw, threads)
    print ('Filtering BAM and removing duplicates on sample {}'.format(sample))
    bam_filtered = sample + '-filtered.bam'
    filter_mapped(bam_raw, bam_filtered, threads)
    bam_dedup = sample + '-dedup.bam'
    remove_duplicates(bam_filtered, bam_dedup, threads)
    bam = sample + '.bam'
    sort(bam_dedup, bam, threads)
    print ('Compute genome coverage for sample {}'.format(sample))
    bam_first_mate = sample + '-mate1.bam'
    first_mate(bam, bam_first_mate, threads)
    bed_first_mate = sample + '-mate1.bed'
    bam_to_bed(bam_first_mate, bed_first_mate)
    count = count_bed(bed_first_mate)
    scale = BASE_SCALE / count
    bed = sample + ".bed"
    bigwig = sample + ".bw"
    genome_coverage(bed_first_mate, bed, sizes, sample, scale)
    bedgraph_to_bigwig(bed, bigwig, sizes)


def download(sample, srr):
    '''Download reads of a sample'''
    fastq = sample + '_1.fastq.gz'
    fastq2 = sample + '_2.fastq.gz'
    srr_output = srr + '_1.fastq.gz'
    srr_output2 = srr + '_2.fastq.gz'
    cmd = ['fastq-dump', '--split-files', '--gzip', srr];
    logging.debug('Running {}'.format(cmd))
    subprocess.call(cmd)
    os.rename(srr_output, fastq)
    if os.path.isfile(srr_output2):
        os.rename(srr_output2, fastq2)


def bwa_index(fasta):
    '''Run BWA on FASTQ files.'''
    bwa_index_cmd = ['bwa', 'index', fasta]
    fasta_indexed = fasta + '.bwt';
    logging.debug('Running {}'.format(bwa_index_cmd))
    subprocess.call(bwa_index_cmd)
    if not os.path.isfile(fasta_indexed):
        raise AssertionError('Error when indexing FASTA ' + fasta)


def bwa(fastq, fastq2, fasta, bam_output, threads=None):
    '''Run BWA on FASTQ files.'''
    sam_output = bam_output + '.sam'
    cmd = ['bwa', 'mem']
    if not threads is None:
        cmd.extend(['-t', str(threads)])
    cmd.extend(['-o', sam_output, fasta, fastq])
    if os.path.isfile(fastq2):
        cmd.append(fastq2)
    logging.debug('Running {}'.format(cmd))
    subprocess.call(cmd)
    if not os.path.isfile(sam_output):
        raise AssertionError('Error when running BWA with command ' + bwa_cmd)
    cmd = ['samtools', 'view', '-b']
    if not threads is None:
        cmd.extend(['--threads', str(threads-1)])
    cmd.extend(['-o', bam_output, sam_output])
    logging.debug('Running {}'.format(cmd))
    subprocess.call(cmd)
    if not os.path.isfile(bam_output):
        raise AssertionError('Error when converting SAM ' + sam_output + ' to BAM ' + bam_output)
    os.remove(sam_output)


def filter_mapped(bam_input, bam_output, threads=None):
    '''Filter BAM file to remove poorly mapped sequences.'''
    cmd = ['samtools', 'view', '-f', '2', '-F', '2048', '-b']
    if not threads is None:
        cmd.extend(['--threads', str(threads-1)])
    cmd.extend(['-o', bam_output, bam_input])
    logging.debug('Running {}'.format(cmd))
    subprocess.call(cmd)
    if not os.path.isfile(bam_output):
        raise AssertionError('Error when filtering BAM ' + bam_input)


def remove_duplicates(bam_input, bam_output, threads=None):
    '''Remove duplicated sequences from BAM file.'''
    fixmate_output = bam_input + '.fix'
    cmd = ['samtools', 'fixmate', '-m']
    if not threads is None:
        cmd.extend(['--threads', str(threads-1)])
    cmd.extend([bam_input, fixmate_output])
    logging.debug('Running {}'.format(cmd))
    subprocess.call(cmd)
    if not os.path.isfile(fixmate_output):
        raise AssertionError('Error when fixing duplicates in BAM ' + bam_input)
    sort_output = bam_input + '.sort'
    cmd = ['samtools', 'sort']
    if not threads is None:
        cmd.extend(['--threads', str(threads-1)])
    cmd.extend(['-o', sort_output, fixmate_output])
    logging.debug('Running {}'.format(cmd))
    subprocess.call(cmd)
    if not os.path.isfile(sort_output):
        raise AssertionError('Error when sorting BAM ' + fixmate_output)
    os.remove(fixmate_output)
    cmd = ['samtools', 'markdup', '-r']
    if not threads is None:
        cmd.extend(['--threads', str(threads-1)])
    cmd.extend([sort_output, bam_output])
    logging.debug('Running {}'.format(cmd))
    subprocess.call(cmd)
    if not os.path.isfile(bam_output):
        raise AssertionError('Error when removing duplicates from BAM ' + bam_output)
    os.remove(sort_output)


def sort(bam_input, bam_output, threads=None):
    '''Sort BAM file.'''
    cmd = ['samtools', 'sort']
    if not threads is None:
        cmd.extend(['--threads', str(threads-1)])
    cmd.extend(['-o', bam_output, bam_input])
    logging.debug('Running {}'.format(cmd))
    subprocess.call(cmd)
    if not os.path.isfile(bam_output):
        raise AssertionError('Error when sorting BAM ' + bam_input)


def first_mate(bam_input, bam_output, threads=None):
    '''Remove second mate from BAM file.'''
    cmd = ['samtools', 'view', '-f', '64', '-b']
    if not threads is None:
        cmd.extend(['--threads', str(threads-1)])
    cmd.extend(['-o', bam_output, bam_input])
    logging.debug('Running {}'.format(cmd))
    subprocess.call(cmd)
    if not os.path.isfile(bam_output):
        raise AssertionError('Error when removing second mate from BAM ' + bam_input)


def bam_to_bed(bam, bed):
    '''Convert BAM file to BED.'''
    cmd = ['bedtools', 'bamtobed', '-i', bam]
    logging.debug('Running {}'.format(cmd))
    with open(bed, "w") as outfile:
        subprocess.call(cmd, stdout=outfile)
    if not os.path.isfile(bed):
        raise AssertionError('Error when converting BAM ' + bam + ' to BED')


def count_bed(bed, strand=None):
    '''Counts number of entry in BED, can be limited to a specific strand.'''
    count = 0
    file = open(bed, 'r')
    for line in file:
        if line.startswith('track') or line.startswith('browser') or line.startswith('#'):
            continue
        if strand is None:
            count += 1
        else:
            columns = line.rstrip('\r\n').split('\t')
            if len(columns) >= 6 and columns[5] == strand:
                count += 1
    file.close()
    return count


def genome_coverage(bed_input, bed_output, sizes, sample, scale=None, strand=None):
    '''Compute genome coverage.'''
    coverage_output = bed_input + '.cov'
    cmd = ['bedtools', 'genomecov', '-bg', '-5', '-i', bed_input, '-g', sizes]
    if not scale is None:
        cmd.extend(['-scale', str(scale)]) 
    if not strand is None:
        cmd.extend(['-strand', strand]) 
    logging.debug('Running {}'.format(cmd))
    with open(coverage_output, "w") as outfile:
        subprocess.call(cmd, stdout=outfile)
    if not os.path.isfile(coverage_output):
        raise AssertionError('Error when computing genome coverage on ' + bed_input)
    sort_output = bed_input + '.sort'
    cmd = ['bedtools', 'sort', '-i', coverage_output]
    logging.debug('Running {}'.format(cmd))
    with open(sort_output, "w") as outfile:
        subprocess.call(cmd, stdout=outfile)
    if not os.path.isfile(sort_output):
        raise AssertionError('Error when sorting BED ' + coverage_output)
    os.remove(coverage_output)
    track = 'track type=bedGraph name="' + sample
    if not strand is None:
        track += ' Minus' if strand == '-' else ' Plus'
    with open(sort_output, "r") as infile:
        with open(bed_output, "w") as outfile:
            outfile.write(track + '\n')
            for line in infile:
                outfile.write(line)
    os.remove(sort_output)


def bedgraph_to_bigwig(bed, bigwig, sizes):
    '''Converts bedgraph file to bigwig.'''
    cmd = ['bedGraphToBigWig', bed, sizes, bigwig]
    logging.debug('Running {}'.format(cmd))
    subprocess.call(cmd)
    if not os.path.isfile(bigwig):
        raise AssertionError('Error when converting BED to BIGWIG ' + bed)


if __name__ == '__main__':
    main()
