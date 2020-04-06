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


rDNA_read_ids = []
id2direction = {}
with open('reads_visualized.txt') as f:
    for item in itertools.zip_longest(*[iter(f)]*2):
        read_id = item[0].split('@')[1].strip()
        rDNA_read_ids.append(read_id)
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
with open('megalodon_cas9/basecalls.fasta') as f:
    for item in itertools.zip_longest(*[iter(f)]*2):
        read_id = item[0][1:].strip()
        read = item[1].strip()
        id2read[read_id] = read

id2scores = {}
with h5py.File('megalodon_cas9/basecalls.modified_base_scores.hdf5') as f:
    for read_id in rDNA_read_ids:
        mc_data = np.array(f['Reads'][read_id])[:,1]
        ma_data = np.array(f['Reads'][read_id])[:,0]
        read = id2read[read_id]
        cpgs = find_cpg(read)
        scores_c = []
        for i in cpgs:
            # each item is [position, methylation_score]
            scores_c.append([i, mc_data[i]])
        id2scores[read_id] = scores_c


def calculate_meth_stats(read_id, start, end):
    scores = [i[1] for i in id2scores[read_id] if start*200 < i[0] < end*200]
    return np.mean(scores)


def analyze_coding(read_id, read, quality, split_length):
    print(read_id)
    allowed_dist = 10000 / split_length
    saminfo = utilities.split_mapping_and_sam_analysis(split_length, read_id, read, quality, 'reference/rDNA_coding_with_prom.fasta')
    codings = []
    coding = []
    prev_n = 0
    for n, line in enumerate(saminfo):
        if line[4] > 70:
            if (n - prev_n) < allowed_dist:
                coding.append(n)
            else:
                if len(coding) > 5:
                    codings.append(coding)
                coding = []
            prev_n = n
    if len(coding) > 5:
        codings.append(coding)
    print(len(codings))
    for coding in codings:
        calculate_meth_stats(read_id, coding[0], coding[-1])


split_length = 200
for read_id in rDNA_read_ids:
    if not(35000 < len(id2read[read_id]) < 55000):
        continue
    read = id2read[read_id]
    if id2direction[read_id] == '-':
        read = str(Seq(read).reverse_complement())
    analyze_coding(read_id, read, 'J'*len(read), split_length)
    continue
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
    if id2direction[read_id] == '+':
        for score in id2scores[read_id]:
            x.append(score[0])
            y.append(math.e**score[1] * 10000)
    if id2direction[read_id] == '-':
        for score in id2scores[read_id]:
            x.append(read_len - score[0])
            y.append(math.e**score[1] * 10000)
    ax.bar(x, y, width=100, color = 'mediumblue')
    ax.add_collection(lc)
    ax.set_yticks((-10000, 0, 10000, 20000, 30000, 40000))
    ax.set_yticklabels(('unmapped', 0, 10000, 20000, 30000, 40000))
    ax.autoscale()
    ax.set_ylim([-12000, 46000])
    plt.savefig('methylation_figs/' + read_id + '.png', dpi=300)
