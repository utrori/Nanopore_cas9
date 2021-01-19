import h5py
import numpy as np
import collections
import matplotlib.pyplot as plt
import glob


class Perread(object):
    def __init__(self, filename):
        self.filename = filename
        with h5py.File(filename) as f:
            self.blocks = list(f['Statistic_Blocks/'].keys())

    def calc_met_stats(self):
        self.get_met_stats()
        met_summary = collections.defaultdict(list)
        for i in self.met_stats:
            met_summary[i[0]].append(i[1])
        for n in range(10000, 24000):
            scores = met_summary.get(n, [0])
            stat = len([i for i in scores if i > -0.5])
            if stat/len(scores) < 0.2:
                print('{}\t{}\t{}'.format(str(n), len(scores), stat))

    def get_met_stats(self):
        with h5py.File(self.filename) as f:
            self.met_stats = []
            for b in self.blocks:
                self.met_stats += self._get_metstat_in_block(f['Statistic_Blocks/' + b])
        return self.met_stats

    def _get_metstat_in_block(self, block):
        temp_met_stats = []
        for stat in block['block_stats']:
            pos, score, id_val = stat
            temp_met_stats.append((pos, score))
        return temp_met_stats

    def get_id2met(self):
        with h5py.File(self.filename) as f:
            self.id2met = {}
            for b in self.blocks:
                self.id2met.update(self._get_id2met_in_block(f['Statistic_Blocks/' + b]))
        return self.id2met

    def _get_id2met_in_block(self, block):
        val2read_id = {}
        id2met_sites = collections.defaultdict(list)
        read_ids = block['read_ids'][()]
        for n, rid in enumerate(read_ids):
            val2read_id[n] = rid.decode()
        scores = []
        poss = []
        for stat in block['block_stats']:
            pos, score, id_val = stat
            id2met_sites[val2read_id[id_val]].append((pos, score))
            scores.append(score)
            poss.append(pos)
        print((max(poss), min(poss)))
        return id2met_sites


def compare_met_stats(stats1, stats2):
    summary1 = collections.defaultdict(list)
    summary2 = collections.defaultdict(list)
    for i in stats1:
        summary1[i[0]].append(i[1])
    for i in stats2:
        summary2[i[0]].append(i[1])
    for n in range(10000, 24000):
        scores = summary1.get(n, [0])
        stat = len([i for i in scores if i > -1])
        if stat/len(scores) < 0.4:
            print('{}\t{}\t{}'.format(str(n+1), len(scores), stat))
            scores2 = summary2[n]
            stat2 = len([i for i in scores2 if i > -1])
            print('{}\t{}\t{}'.format(str(n+1), len(scores2), stat2))



class Resquiggled_f5(object):
    def _find_basecall_with_mod(self, h5file):
        for i in range(3):
            bc_template = h5file['Analyses/Basecall_1D_00' + str(i) + '/BaseCalled_template']
            if 'ModBaseProbs' in bc_template.keys():
                return str(i)

    def __init__(self, filename):
        self.filename = filename
        with h5py.File(filename) as f:
            read_name = list(f['Raw/Reads'].keys())[0]
            self.read_id = f['Raw/Reads/' + read_name].attrs['read_id'].decode()
            mod_loc = self._find_basecall_with_mod(f)
            fastq = f['/Analyses/Basecall_1D_00' + mod_loc + '/BaseCalled_template/Fastq']\
                    [()].decode().strip().split('\n')
            self.seq = fastq[1]
            self.res_bases = ''
            events = f['/Analyses/RawGenomeCorrected_000/BaseCalled_template/Events'][()]
            for e in events:
                self.res_bases += e[4].decode()


if __name__ == '__main__':
    r = Resquiggled_f5('/home/yutaro/nanopore/cas9/1015_head_met_basecalled/workspace/0b2e43ec-e121-4611-a75a-61d66317eee3.fast5')
    p = Perread('1015_met_cpg_perbase.CpG.tombo.per_read_stats')
    p2 = Perread('1015_unmet_cpg_perbase.CpG.tombo.per_read_stats')
    met_stat1 = p.get_met_stats()
    met_stat2 = p2.get_met_stats()
    compare_met_stats(met_stat1, met_stat2)
    quit()
    id2met = p.get_id2met()
    themet = id2met[r.read_id]
    poss = [i[0] for i in themet]
    for pos in poss:
        print('---------------')
        print(r.seq[pos])
        print(r.res_bases[pos])
