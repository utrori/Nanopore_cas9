import utilities
import subprocess
import os
import sys
import search_rDNA_reads
import glob
import itertools
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd


def unmapped_coords(coords):
    unmapped_dist = []
    flag = 0
    true_coords = []
    for coord in coords:
        if flag == 0:
            if coord != 0:
                flag = 1
        if flag == 1:
            true_coords.append(coord)
    current = true_coords[0]
    for coord in true_coords:
        if coord == 0:
            unmapped_dist.append(current)
        else:
            current = coord
    return unmapped_dist


def mapping():
    files = glob.glob('/home/yutaro/nanopore/clive/clive_rDNA_reads.fastq')
    FNULL = open(os.devnull, 'w')
    split_length = 200
    unmapped_dist = []
    count = 0

    for fastq_file in files:
        length_cutoff = 35000
        with open(fastq_file) as f:
            for data in itertools.zip_longest(*[iter(f)]*4):
                count += 1
                if count % 100 == 0:
                    print(count)
                header = data[0].strip()
                read = data[1].strip()
                if len(read) < length_cutoff:
                    continue
                quality = data[3].strip()
                utilities.make_temp_fastq(split_length, header, read, quality)
                subprocess.run('bwa mem -M -x ont2d -t 5 '
                               '/home/yutaro/nanopore/clive/'
                               'rDNA_index/humRibosomal.fa temp_files/temp_fastq.fastq > '
                               'temp_files/temp_sam.sam', shell=True, stdout=FNULL,
                               stderr=subprocess.STDOUT)
                with open('temp_files/temp_sam.sam') as samf:
                    samdata = samf.readlines()[2:]
                    ASs = []
                    coords = []
                    for line in samdata:
                        row = line.split()
                        AS = line.split('AS:i:')[1].split()[0]
                        ASs.append(int(AS))
                        coord = row[3]
                        coords.append(int(coord))
                    if np.mean(ASs) > 70:
                        unmapped_dist.extend(unmapped_coords(coords))
    pd.to_pickle(unmapped_dist, 'unmapped.pkl')
    return unmapped_dist

if __name__ == '__main__':
    #unmapped_dist = mapping()
    unmapped_dist = pd.read_pickle('unmapped.pkl')
    plt.xticks([0, 1000,2000,3000,4000,5000,6000,7000,8000, 9000, 10000])
    plt.hist(unmapped_dist, bins=500)
    plt.savefig('unmapped/clive_duplicate2.png', dpi=300)
