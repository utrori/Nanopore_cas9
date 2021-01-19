import h5py
import shutil
import glob
import os
import utilities
import fast5class
import numpy as np

with open('settings.txt') as f:
    ref = f.readline().split()[1].strip()

def search_rDNA_reads(in_dir_name, rDNA_file_name):
    if in_dir_name[-1] == '/':
        in_dir_name = in_dir_name[:-1]
    fast5files = glob.glob(in_dir_name + '/0/*.fast5')
    string = ''
    for n, fast5 in enumerate(fast5files):
        if n % 1000 == 0:
            print(n)
        read = fast5class.Read(fast5, ref)
        if 4000 < read.length:
            scores = read.sam_summary[:,4].astype(np.uint16)
            if len(np.nonzero(scores)[0]) / len(scores) > 0.8:
                string += read.read_id + '\n'
    with open(rDNA_file_name, 'w') as fw:
        fw.write(string)


def move_rDNA_reads():
    fast5_ids = []
    with open('1020_rDNA_reads.txt') as f:
        for line in f:
            fast5_ids.append(line.strip())
    new_dir_name = '1020_rDNA_singles'
    if os.path.exists(new_dir_name):
        shutil.rmtree(new_dir_name)
    os.mkdir(new_dir_name)
    fast5files = glob.glob('1020_single_fast5s/*/*.fast5')
    for fast5 in fast5files:
        if any((i in fast5 for i in fast5_ids)):
            shutil.copy(fast5, new_dir_name)


if __name__ == '__main__':
    search_rDNA_reads('1215_BSL2KA_bc_fast5s/', '1215_BSL_rDNAs.txt')
