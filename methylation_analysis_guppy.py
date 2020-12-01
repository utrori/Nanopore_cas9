# coding: utf-8
import glob
import math
import collections
from visualize_nanopore_read import plot_read_structure
from ont_fast5_api.fast5_interface import get_fast5_file
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


rDNA_read_ids = []
with open('rDNA_containing_reads.fastq') as f:
    for item in itertools.zip_longest(*[iter(f)]*4):
        read_id = item[0].split()[0].strip()[1:]
        rDNA_read_ids.append(read_id)
        summary = 0
        read_vis = item[1]

def guppy_fast5_extraction(filename):
    data_set = []
    with get_fast5_file(filename, mode='r') as f:
        for read_id in f.get_read_ids():
            read = f.get_read(read_id)
            latest_basecall = read.get_latest_analysis('Basecall_1D')
            fastq = read.get_analysis_dataset(latest_basecall, 
            'BaseCalled_template/Fastq')
            mod_base_table = read.get_analysis_dataset(latest_basecall, 
            'BaseCalled_template/ModBaseProbs')
            data_set.append((read_id, fastq, mod_base_table))
    return data_set


def make_methylation_summary(read, mod_base_table):
    cpgs = find_cpg(read)
    mod_scores = []
    for cpg in cpgs:
        scores = mod_base_table[cpg]
        mod_scores.append((cpg, scores[3]/255))
    return mod_scores


def calculate_meth_stats(mod_scores):
    meth_stats = collections.defaultdict(list)
    bin_met= 200
    for score in mod_scores:
        meth_stats[score[0]//200].append(score[1])
    ret_met = {}
    for n, score in meth_stats.items():
        ret_met[n] = np.mean(score)
    return ret_met


dirpath = '293_met/workspace/'

contexts = []
split_length = 200
for fast5 in glob.glob(dirpath + '*.fast5'):
    dataset = guppy_fast5_extraction(fast5)
    for data in dataset:
        read_id = data[0]
        if read_id not in rDNA_read_ids:
            continue
        split_fastq = data[1].split()
        read = split_fastq[6]
        quality = split_fastq[8]
        mod_base_table = data[2]
        cpg_scores = make_methylation_summary(read, mod_base_table)
        averaged_cpg = calculate_meth_stats(cpg_scores)
        utilities.split_mapping_and_sam_analysis(split_length, read_id, read, quality, '../clive/rDNA_index/humRibosomal.fa')
        lc = plot_read_structure('test', split_length, 0)
        fig = plt.figure()
        plt.subplots_adjust(left=0.2)
        ax = fig.add_subplot()
        x = []
        y = []
        for n, score in averaged_cpg.items():
            x.append(n * 200)
            y.append(score * 10000)
        xa = []
        ya = []
        for n, base in enumerate(mod_base_table):
            if base[1] != 0 and read[n] == 'A':
                contexts.append(read[n-1:n+3])
                xa.append(n)
                ya.append(base[1] * 40)
        ax.bar(x, y, width=200, color = 'mediumblue', zorder=0)
        ax.bar(xa, ya, width=200, color = 'red', zorder=0)
        ax.add_collection(lc)
        ax.set_yticks((-10000, 0, 10000, 20000, 30000, 40000))
        ax.set_yticklabels(('unmapped', 0, 10000, 20000, 30000, 40000))
        ax.autoscale()
        ax.set_ylim([-12000, 46000])
        plt.savefig('293_visualize/' + read_id + '.png', dpi=300)
        plt.close()
with open('6mA_context_293.txt', 'w') as fw:
    for i in contexts:
        fw.write(i + '\n')
