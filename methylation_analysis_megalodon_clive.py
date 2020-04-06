# coding: utf-8
import math
from visualize_nanopore_read import plot_read_structure
import os
import search_rDNA_reads
from rDNA_structure_of_long_reads import easy_flag
import subprocess
from Bio.Seq import Seq
import h5py
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import itertools
import utilities


def find_cpg(read):
    cpg_sites = []
    read = read.upper()
    for n in range(len(read)-1):
        if read[n] == 'C':
            if read[n+1] == 'G':
                cpg_sites.append(n)
    return cpg_sites


filenames = os.listdir('clive/')
rDNA_read_ids = [i.split('.')[0] for i in filenames if 'fast5' in i]

id2direction = {}
with open('../clive/reads_visualized.txt') as f:
    for item in itertools.zip_longest(*[iter(f)]*2):
        read_id = item[0].split('@')[1].strip()
        summary = 0
        read_vis = item[1]
        for item in read_vis:
            if item =='￫':
                summary += 1
            elif item =='￩':
                summary -= 1
        if summary > 0:
            id2direction[read_id] = '+'
        elif summary < 0:
            id2direction[read_id] = '-'

id2read = {}
with open('megalodon_results/basecalls.fasta') as f:
    for item in itertools.zip_longest(*[iter(f)]*2):
        read_id = item[0][1:].strip()
        read = item[1].strip()
        id2read[read_id] = read

id2scores = {}
with h5py.File('megalodon_results/basecalls.modified_base_scores.hdf5') as f:
    for read_id in rDNA_read_ids:
        mc_data = np.array(f['Reads'][read_id])[:,1]
        ma_data = np.array(f['Reads'][read_id])[:,0]
        read = id2read[read_id]
        cpgs = find_cpg(read)
        scores_c = []
        for i in cpgs:
            scores_c.append([i, mc_data[i]])
        id2scores[read_id] = scores_c

split_length = 200
for read_id in rDNA_read_ids:
    read = id2read[read_id]
    search_rDNA_reads.make_temp_fastq(split_length, read_id, read, 'J'*len(read))
    FNULL = open(os.devnull, 'w')
    subprocess.run('bwa mem -M -x ont2d -t 5 reference/rDNA_for_cas9_no_margin.fasta '
                               'temp_files/temp_fastq.fastq > temp_files/'
                               'single_split_mapped.sam',
                               shell=True, stdout=FNULL,
                               stderr=subprocess.STDOUT)
    lc = plot_read_structure('test', split_length, 9400)
    read_len = len(read)
    fig = plt.figure()
    plt.subplots_adjust(left=0.2)
    ax = fig.add_subplot()
    x = []
    y = []
    for score in id2scores[read_id]:
        x.append(score[0])
        y.append(math.e**score[1] * 10000)
    ax.bar(x, y, width=100, color = 'mediumblue')
    ax.add_collection(lc)
    ax.set_yticks((-10000, 0, 10000, 20000, 30000, 40000))
    ax.set_yticklabels(('unmapped', 0, 10000, 20000, 30000, 40000))
    ax.autoscale()
    ax.set_ylim([-12000, 46000])
    plt.savefig('methylation_figs_clive/' + read_id + '.png', dpi=300)
    plt.close()
