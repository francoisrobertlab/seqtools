from distutils.command.check import check
import logging
import os
import subprocess

import click
from seqtools.txt import Parser

genes1 = [columns[3] for columns in Parser.columns('top-q1.txt')]
genes2 = [columns[3] for columns in Parser.columns('top-q2.txt')]
genes3 = [columns[3] for columns in Parser.columns('top-q3.txt')]
genes4 = [columns[3] for columns in Parser.columns('top-q4.txt')]
genes5 = [columns[3] for columns in Parser.columns('top-q5.txt')]
genestop10 = [columns[3] for columns in Parser.columns('Top10percentTxbdGenes.txt')]
allgenes = Parser.first('sgdGeneAndOthers_UTR_TIF-seq_sacCer3_july_2018.txt')
for gene in genes1:
    if gene in genes2:
        print ('gene {} in q1 and q2'.format(gene))
    if gene in genes3:
        print ('gene {} in q1 and q3'.format(gene))
    if gene in genes4:
        print ('gene {} in q1 and q4'.format(gene))
    if gene in genes5:
        print ('gene {} in q1 and q5'.format(gene))
for gene in genes2:
    if gene in genes3:
        print ('gene {} in q1 and q3'.format(gene))
    if gene in genes4:
        print ('gene {} in q1 and q4'.format(gene))
    if gene in genes5:
        print ('gene {} in q1 and q5'.format(gene))
for gene in genes3:
    if gene in genes4:
        print ('gene {} in q1 and q4'.format(gene))
    if gene in genes5:
        print ('gene {} in q1 and q5'.format(gene))
for gene in genes4:
    if gene in genes5:
        print ('gene {} in q1 and q5'.format(gene))
for gene in genestop10:
    if gene not in genes1:
        print ('gene {} in top10 but not in q1'.format(gene))
    if gene in genes2:
        print ('gene {} in top10 and q2'.format(gene))
    if gene in genes3:
        print ('gene {} in top10 and q3'.format(gene))
    if gene in genes4:
        print ('gene {} in top10 and q4'.format(gene))
    if gene in genes5:
        print ('gene {} in top10 and q5'.format(gene))
for gene in allgenes:
    if gene not in genes1 and gene not in genes2 and gene not in genes3 and gene not in genes4 and gene not in genes5:
        print ('gene {} not in any quantiles'.format(gene))
