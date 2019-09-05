import glob
import logging
import os
import shutil
import subprocess

import FullAnalysis
import SplitBed
import click


@click.command()
@click.option('--samples', '-s', type=click.Path(exists=True), default='samples.txt',
              help='Sample names listed one sample name by line.')
@click.option('--parameters', '-p', type=click.Path(exists=True), default='parameters.txt',
              help='VAP parameters file.')
def main(samples, parameters):
    '''Run VAP on samples.'''
    logging.basicConfig(filename='debug.log', level=logging.DEBUG, format='%(asctime)s %(levelname)-8s %(message)s', datefmt='%Y-%m-%d %H:%M:%S')
    samples_names = FullAnalysis.first_column(samples)
    for sample in samples_names:
        vap(sample, parameters)


def vap(sample, parameters):
    '''Run VAP on a single sample.'''
    print ('Running VAP on sample {}'.format(sample))
    output = sample + '-vap-output';
    if not os.path.exists(output):
        os.mkdir(output);
    sample_parameters = output + '/parameters.txt';
    splits = SplitBed.splits(sample)
    beds = [split + '.bed' for split in splits]
    genes = parse_genes(parameters)
    create_parameters(beds, output, parameters, sample_parameters);
    cmd = ['./vap']
    if os.name == 'nt':
        cmd = ['./vap.exe']
    cmd.extend(['-p', sample_parameters])
    logging.debug('Running {}'.format(cmd))
    subprocess.call(cmd)
    try:
        splits_values = parse_heatmap_values(splits, output)
    except AssertionError as exception:
        raise AssertionError('Error when running VAP with parameters ' + sample_parameters, exception)
    merged_heatmap = sample + '-heatmap.txt';
    create_heatmap(sample, genes, splits, splits_values, merged_heatmap)
    shutil.rmtree(output)


def create_parameters(datasets, output_folder, parameters_input, parameters_output):
    with open(parameters_input, 'r') as infile:
        with open(parameters_output, 'w') as outfile:
            for parameters_line in infile:
                if parameters_line.startswith('~~@dataset_path'):
                    for dataset in datasets:
                        dataset_line = '~~@dataset_path=R1:=:' + dataset
                        outfile.write(dataset_line)
                        outfile.write('\n')
                elif parameters_line.startswith('~~@output_directory'):
                    parameters_line = '~~@output_directory=' + output_folder
                    outfile.write(parameters_line)
                    outfile.write('\n')
                elif parameters_line.startswith('~~@prefix_filename='):
                    index = parameters_line.index('=')
                    outfile.write(parameters_line[:index + 1])
                    outfile.write('\n')
                else:
                    outfile.write(parameters_line)


def parse_genes(parameters):
    with open(parameters, 'r') as infile:
        for parameters_line in infile:
            if parameters_line.startswith('~~@refgroup_path='):
                index = parameters_line.index(':=:') + 3
                genes_file = parameters_line.rstrip('\r\n')[index:]
    with open(genes_file, 'r') as infile:
        genes = [line.rstrip('\r\n') for line in infile]
        genes = [line for line in genes if not line.startswith('#')]
    return genes


def parse_heatmap_values(sample_splits, output_folder):
    error_count = 0
    splits_values = {}
    for split in sample_splits:
        splits_values[split] = {}
        split_glob = output_folder + '/ind_data_' + split + '*.txt';
        split_data = glob.glob(split_glob)
        if not split_data:
            error_count += 1
            logging.warning('Cannot open VAP output file {}'.format(split_glob))
            continue
        with open(split_data[0], 'r') as infile:
            index = None
            for line in infile:
                if not line.startswith('#'):
                    columns = line.rstrip('\r\n').split('\t');
                    if index is None:
                        index = columns.index('W0_0')
                    else:
                        gene = columns[0]
                        value = columns[index]
                        splits_values[split][gene] = value
    if error_count == len(sample_splits):
        raise AssertionError('Cannot open any VAP output file')
    return splits_values


def create_heatmap(sample, genes, splits, splits_values, output):
    with open(output, 'w') as outfile:
        outfile.write('UNIQID')
        outfile.write('\t')
        outfile.write('Name')
        for split in splits:
            outfile.write('\t')
            outfile.write(split[len(sample) + 1:])
        outfile.write('\n')
        outfile.write('EWEIGHT')
        outfile.write('\t')
        for split in splits:
            outfile.write('\t')
            outfile.write('1')
        outfile.write('\n')
        for gene in genes:
            outfile.write(gene)
            outfile.write('\t')
            outfile.write(gene)
            for split in splits:
                outfile.write('\t')
                if splits_values[split] and splits_values[split][gene]:
                    outfile.write(splits_values[split][gene])
                else:
                    outfile.write('0')
            outfile.write('\n')


if __name__ == '__main__':
    main()
