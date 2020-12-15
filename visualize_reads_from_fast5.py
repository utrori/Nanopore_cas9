import h5py
import matplotlib.pyplot as plt
import utilities
import numpy as np
import matplotlib.cm as cm
import glob
import os
import shutil

class Read(object):
    def _find_basecall_with_mod(self, h5file):
        for i in range(3):
            bc_template = h5file['Analyses/Basecall_1D_00' + str(i) + '/BaseCalled_template']
            if 'ModBaseProbs' in bc_template.keys():
                return str(i)

    def __init__(self, filename):
        with h5py.File(filename) as f:
            read_name = list(f['Raw/Reads'].keys())[0]
            self.read_id = f['Raw/Reads/' + read_name].attrs['read_id'].decode()
            self.raw_duration = f['Raw/Reads/' + read_name].attrs['duration']
            self.raw = f['Raw/Reads/'+ read_name + '/Signal'][:]
            mod_loc = self._find_basecall_with_mod(f)
            fastq = f['/Analyses/Basecall_1D_00' + mod_loc + '/BaseCalled_template/Fastq']\
                    [()].decode().strip().split('\n')
            self.seq = fastq[1]
            self.quality = fastq[3]
            #fastq length! not raw signal length
            self.length = len(self.seq)
            self.guppy_mbt = \
                    f['/Analyses/Basecall_1D_00' + mod_loc + '/BaseCalled_template/ModBaseProbs'][:]
            self.split_length = 300
            sam_summary = np.array(utilities.split_mapping_and_sam_analysis(self.split_length, 'temp', 
                    self.seq, self.quality, ref))

    def _cpg_methylation_average(self):
        x = []
        y = []
        window = 100
        cpg_met_scores = self.guppy_mbt[:,3]
        for n, item in enumerate([cpg_met_scores[i:i+window] for i in range(0, len(cpg_met_scores), window)]):
            x.append(n * window)
            y.append(len([i for i in item if i>180]) / window * 70000)
        return x, y

    def plot_structure(self, ref, savedir):
        utilities.split_mapping_and_sam_analysis(self.split_length, self.read_id, self.seq, self.quality, ref)
        lc = utilities.plot_read_structure(self.read_id, self.split_length)
        fig = plt.figure()
        plt.subplots_adjust(left=0.2)
        ax = fig.add_subplot()
        x, y = self._cpg_methylation_average()
        """
        for n, item in enumerate(self.guppy_mbt[:,1]):
            if item > 120:
                ax.bar(n, item * 10000/255, width=100, color='red', zorder=3)
        """
        ax.bar(x, y, width=100, color = 'mediumblue', zorder=0)
        ax.add_collection(lc)
        ax.set_yticks((-10000, 0, 10000, 20000, 30000, 40000))
        ax.set_yticklabels(('unmapped', 0, 10000, 20000, 30000, 40000))
        ax.autoscale()
        ax.set_ylim([-12000, 46000])
        if savedir[-1] == '/':
            savedir = savedir[:-1]
        plt.savefig(savedir + '/' + self.read_id + '.png', dpi=300)
        plt.close()


if __name__ == '__main__':
    ref = 'rDNA_index/humRibosomal.fa'
    fast5_dir = '1215_BSL2KA_bc_fast5s/0'
    fast5files = glob.glob(fast5_dir + '/*.fast5')
    print(fast5files)
    reads = []
    savedir = '1215_plot'
    if os.path.exists(savedir):
        shutil.rmtree(savedir)
    os.mkdir(savedir)
    rDNA_ids = []
    with open('1215_BSL_rDNAs.txt') as f:
        for line in f:
            rDNA_ids.append(line.strip())
    for fast5 in fast5files:
        if fast5.split('.')[0].split('/')[-1] in rDNA_ids:
            r = Read(fast5)
            r.plot_structure(ref, savedir)
