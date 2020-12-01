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

def mapping_find_end():
    files = glob.glob('/var/lib/minknow/data/201015_cas/517/20201015_0457_MN32877_FAO42372_5d01cfab/fastq_pass/*.fastq')
    FNULL = open(os.devnull, 'w')
    split_length = 200
    end_dist = []

    for fastq_file in files:
        length_cutoff = 1000
        with open(fastq_file) as f:
            for data in itertools.zip_longest(*[iter(f)]*4):
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
                    if np.mean(ASs) > 100:
                        if coords[-1] != 0:
                            end_dist.append(coords[-1])
                        elif coords[-2] != 0:
                            end_dist.append(coords[-2])

    pd.to_pickle(end_dist, 'end_dist.pkl')

if __name__ == '__main__':
    end_dist = pd.read_pickle('end_dist.pkl')
    plt.hist(end_dist, bins=50)
    plt.savefig('end_dist/517.png', dpi=300)
