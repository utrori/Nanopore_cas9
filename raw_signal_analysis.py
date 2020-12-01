import matplotlib.pyplot as plt
import numpy as np
from ont_fast5_api.fast5_interface import get_fast5_file
import h5py
import os
import itertools
import glob


rDNA_ids = []
with open('1015_rDNA_reads.fastq') as f:
    for items in itertools.zip_longest(*[iter(f)]*4):
        rDNA_ids.append(items[0].split()[0][1:])


def raw_signal_fast5(filename):
    count = 0
#    with get_fast5_file(filename, 'r') as f:
#        for read_id in f.get_read_ids():
#            if read_id in rDNA_ids:
#                count += 1
#                read = f.get_read(read_id)
#                raw = read.get_raw_data()
#                break
    with h5py.File('/home/yutaro/nanopore/cas9/1015_single_fast5s/0/0a0b9e3a-0731-4087-aa06-819f5b5f5e2a.fast5', 'r') as f:
        raw = f['Raw']
        with h5py.File('temp_fast5s/test_fast5.fast5', 'w') as outf:
            f.copy('', outf)
            outf.flush()


if __name__ == '__main__':
    raw = raw_signal_fast5('/var/lib/minknow/data/201015_cas/517/20201015_0457_MN32877_FAO42372_5d01cfab/fast5_pass/FAO42372_pass_0e8e9c19_0.fast5')
    quit()
    x = [i for i in range(len(raw))]
    plt.plot(x[5000:5700], raw[5000:5700], zorder=0)
    plt.scatter(x[5000:5700], raw[5000:5700], s=1, color='red',zorder=1)
    plt.savefig('temp_files/temp.png', dpi=300)
