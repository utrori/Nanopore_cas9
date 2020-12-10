import h5py
import shutil
import glob
import os
import utilities
import numpy as np

ref = 'rDNA_index/humRibosomal.fa'

def search_rDNA_reads(in_dir_name, rDNA_file_name):
    fast5files = glob.glob(in_dir_name + '/workspace/*.fast5')
    string = ''
    for n, fast5 in enumerate(fast5files):
        if n % 1000 == 0:
            print(n)
        with h5py.File(fast5) as f:
            read_name = list(f['Raw/Reads'].keys())[0]
            read_id = f['Raw/Reads/' + read_name].attrs['read_id'].decode()
            fastq = f['/Analyses/Basecall_1D_000/BaseCalled_template/Fastq']\
                        [()].decode().strip().split('\n')
            seq = fastq[1]
            quality = fastq[3]
            length = len(seq)
            sam_summary = np.array(utilities.split_mapping_and_sam_analysis(300, 'temp', 
                    seq, quality, ref))
            if 40000 < length < 44000:
                scores = sam_summary[:,4].astype(np.uint16)
                if len(np.nonzero(scores)[0]) / len(scores) > 0.8:
                    string += read_id + '\n'
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
    search_rDNA_reads()
