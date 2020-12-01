import numpy as np
import matplotlib.pyplot as plt
import h5py
import itertools
import make_truncated_fast5

rDNA_fastq_file = '/home/yutaro/nanopore/cas9/1015_rDNA_reads.fastq'
id2seq = make_truncated_fast5.rDNA_read_ids2seq(rDNA_fastq_file)

def make_dir_met_summary():
    with open('1015_dir_met_summary.txt', 'w') as fw:
        for rid in id2seq:
            with h5py.File('1015_mod_fast5s/workspace/' + rid + '.fast5') as f:
                fastq = f['Analyses/Basecall_1D_001/BaseCalled_template/Fastq'][()]
                modbasetable = f['Analyses/Basecall_1D_001/BaseCalled_template/ModBaseProbs']
                seq = id2seq[rid]
                length = len(seq)
                if 39000 < length < 44000:
                    direction = make_truncated_fast5.get_direction(seq)
                    if direction == '+':
                        reg1 = [0, 3500]
                        reg2 = [length - 9000, length]
                    elif direction == '-':
                        reg1 = [0, 9000]
                        reg2 = [length - 3500, length]
                    reg1_5mc_ave = np.mean(modbasetable[reg1[0]:reg1[1], 3])
                    reg2_5mc_ave = np.mean(modbasetable[reg2[0]:reg2[1], 3])
                    fw.write(rid + '\t' + direction + '\t' + str(reg1_5mc_ave) + '\t' + str(reg2_5mc_ave) + '\n')

mc_aves = []
with open('1015_dir_met_summary.txt') as f:
    for line in f:
        row = line.strip().split()
        mc_aves.append(float(row[2]))
        mc_aves.append(float(row[3]))

plt.hist(mc_aves, bins=50)
plt.show()
