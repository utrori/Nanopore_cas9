import h5py
import time
import re
import os
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import glob
import utilities
import subprocess
import shutil
from Bio.Seq import Seq
from cigar import Cigar

FNULL = open(os.devnull, 'w')

standard_ref_seq = ''
with open('rDNA_index/humRibosomal.fa') as f:
    for l in f.readlines()[1:]:
        standard_ref_seq += l.strip()


def get_ref_offset(ref):
    with open('temp_files/temp_index.fa', 'w') as fw:
        fw.write('>temp_index\n' + standard_ref_seq[0:200])
    subprocess.run('bwa mem -M -x ont2d -t 7 ' + ref + 
                   ' temp_files/temp_index.fa > temp_files/'
                   'temp.sam',
                   shell=True, stdout=FNULL,
                   stderr=subprocess.STDOUT)
    with open('temp_files/temp.sam') as f:
        return int(f.readlines()[2].split()[3]) - 1


def calc_qual(qual_str):
    qual_vals = []
    for s in qual_str:
        qual_vals.append(ord(s))
    qual_actual_prob = [10**((-i+33)/10) for i in qual_vals]
    smalls = [i for i in qual_vals if i < 37]
    return np.median(qual_actual_prob)


class Read(object):
    def _is_this_healthy_rDNA(self):
        """Check if the read is completely inside rDNA or at the end.
        """
        if self.length < 20000:
            return 0
        mapping_state = []
        for item in self.sam_summary:
            if item[1] != '0':
                mapping_state.append(1)
            else:
                mapping_state.append(0)
        threshold = 0.8
        if sum(mapping_state)/len(mapping_state) > threshold:
            return 1
        else:
            for i in range(1, len(mapping_state) - 50):
                if sum(mapping_state[i:])/len(mapping_state[i:]) > threshold or \
                   sum(mapping_state[:-i])/len(mapping_state[:-i]) > threshold:
                    healthy = 2
        return 0
 
    def _get_direction(self, direction_array):
        minus = 0.1
        plus = 0.1
        for i in direction_array:
            if i == '+':
                plus += 1
            if i == '-':
                minus += 1
        if plus/minus > 0.9:
            return '+'
        elif minus/plus > 0.9:
            return '-'
        else:
            return 0

    def _find_basecall_with_mod(self, h5file):
        for i in range(3):
            bc_template = h5file['Analyses/Basecall_1D_00' + str(i) + '/BaseCalled_template']
            if 'ModBaseProbs' in bc_template.keys():
                return str(i)

    def __init__(self, fast5file, ref):
        #Use fast5s that are basecalled with modification!
        self.fast5path = fast5file
        with h5py.File(fast5file) as f:
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
            self.split_length = 300
            self.length = len(self.seq)
            self.guppy_mbt = \
                    f['/Analyses/Basecall_1D_00' + mod_loc + '/BaseCalled_template/ModBaseProbs'][:]
            self.sam_summary = np.array(utilities.split_mapping_and_sam_analysis(self.split_length, 'temp', 
                    self.seq, self.quality, ref))
            self.direction = self._get_direction(self.sam_summary[:,1])
            self.long_side_len = 9500
            self.short_side_len = 3300

    def is_this_rDNA_read(self):
        n = 0
        for item in self.sam_summary:
            if item[3] != '*':
                n += 1
        if n / len(self.sam_summary) > 0.5:
            return 1
        else:
            return 0

    def get_methylated_bases(self, threshold):
        cpg = self.guppy_mbt[:,3]
        dam = self.guppy_mbt[:,1]
        cpg_pos = np.where(cpg > threshold)
        dam_pos = np.where(dam > threshold)
        """
        Since the outputs of np.where is tuple we need to slice by [0].
        The reason the output is tuple is so that it can work on higher
        dimensional arrays.
        """
        return cpg_pos[0], dam_pos[0]

    def get_coding_methylation_pros(self, met_type):
        """
        get coding methylation propotion for visualization
        """
        cpg_pos, dam_pos = self.get_methylated_bases(180)
        if met_type.lower()== 'cpg':
            met_pos = cpg_pos
        if met_type.lower() == 'dam':
            met_pos = dam_pos
        short_side = self.short_side_len
        long_side = self.long_side_len
        if self.direction == '+':
            head = [i for i in met_pos if 0 < i < short_side]
            tail = [i for i in met_pos if self.length-long_side < i < self.length]
            return len(head)/short_side, len(tail)/long_side
        elif self.direction == '-':
            head = [i for i in met_pos if 0 < i < long_side]
            tail = [i for i in met_pos if self.length-short_side < i < self.length]
            return len(head)/long_side, len(tail)/short_side

    def get_unmapped_section(self):
        print(self.read_id)
        print(calc_qual(self.quality))
        pos = self.sam_summary[:,2].astype(np.uint32)
        for n, p in enumerate(pos):
            if p == 0:
                um_qual = self.quality[n*self.split_length:(n+1)*self.split_length]
                if um_qual:
                    print(calc_qual(um_qual))

    def __coding_coordinate(self):
        """
        get coding coordinate of read, not reference
        """
        region1 = self.long_side_len
        region2 = self.short_side_len
        length = len(self.seq)
        if self.direction == '+':
            a_s = 0
            a_e = region2
            b_s = self.length - region1
            b_e = self.length - 1
        elif self.direction == '-':
            a_s = 0
            a_e = region1
            b_s = self.length - region2
            b_e = self.length - 1
        return (a_s, a_e, b_s, b_e)

    def _truncate_fast5_file(self, outname, start, end, read_id=''):
        filename = self.fast5path
        length = end - start
        with h5py.File(filename, 'r') as f:
            with h5py.File(outname, 'w') as fw:
                f.copy('Raw', fw)
                f.copy('UniqueGlobalKey', fw)
                read_name = list(f['Raw/Reads'].keys())[0]
                if read_id:
                    fw['Raw/Reads/'+ read_name].attrs.modify('read_id', read_id.encode())
                fw['/'].attrs.create('file_version', 2, dtype='double')
                fw['Raw/Reads/' + read_name].attrs['duration'] = length
                signal = f['Raw/Reads/'+ read_name + '/Signal']
                del fw['Raw/Reads/'+ read_name + '/Signal']
                temp = fw['Raw/Reads/' + read_name].create_dataset('Signal', \
                        (length,), dtype='<i2', maxshape=(None,), \
                        chunks=(length,), compression=1)
                temp[:] = signal[start:end]
                fw.flush()

    def _get_raw_pos(self, coordinate):
        signal_length = len(self.raw)
        return int(signal_length/self.length * coordinate)

    def get_position_in_ref(self, pos, ref):
        around_seq = self.seq[pos-50:pos+50]
        around_qual = self.quality[pos-50:pos+50]
        with open('temp_files/temp_fastq.fastq', 'w') as fw:
            fw.write('@temp\n' + around_seq + '\n+\n' + 
                     around_qual)
        subprocess.run('bwa mem -M -x ont2d -t 7 ' + ref + 
                       ' temp_files/temp_fastq.fastq > temp_files/'
                       'temp.sam',
                       shell=True, stdout=FNULL,
                       stderr=subprocess.STDOUT)
        with open('temp_files/temp.sam') as f:
            row = f.readlines()[2].strip().split()
            if row[2] == '*':
                return None
            cigar = row[5]
            ref_pos = int(row[3])
            c = Cigar(cigar)
            split_cigar = ''
            for i in c.items():
                split_cigar += i[0] * i[1]
            shift = 0
            current = 0
            for l in split_cigar:
                current += 1
                if l == 'I':
                    shift -= 1
                elif l == 'D':
                    shift += 1
                    current -= 1
                if current == 50:
                    ref_coordinate = ref_pos + 49 + shift
                    break
            return ref_coordinate

    def truncate_and_separate_by_dam(self, start, end, fast5_dir, control_dir, offset):
        #separate reads by the dam status within (start, end) reference coordinate
        coding_coords = self.__coding_coordinate()
        if end < 9200 or start > 20000:
            if self.direction == '+':
                target_side = 'tail'
            else:
                target_side = 'head'
        elif start > 10000:
            if self.direction == '+':
                target_side = 'head'
            else:
                target_side = 'tail'
        else:
            return None
        raw_coords = [self._get_raw_pos(i) for i in coding_coords]
        if target_side == 'head':
            raw_pos = (raw_coords[0], raw_coords[1])
        else:
            raw_pos = (raw_coords[2], raw_coords[3])
        dam_met_pos = self.get_methylated_bases(180)[1]
        ref_dam_pos = [self.get_position_in_ref(i, ref) for i in dam_met_pos]
        ref_dam_pos = [i for i in ref_dam_pos if i != None]
        if not ref_dam_pos:
            self._truncate_fast5_file(control_dir + '/' + self.read_id + '.fast5', 
                                      raw_pos[0], raw_pos[1])
        elif any(((start + offset) % 42999 < i < (end + offset) % 42999 for i in ref_dam_pos)):
            self._truncate_fast5_file(fast5_dir + '/' + self.read_id + '.fast5', 
                                      raw_pos[0], raw_pos[1])
        else:
            self._truncate_fast5_file(control_dir + '/' + self.read_id + '.fast5', 
                                      raw_pos[0], raw_pos[1])

    def truncate_and_separation_by_cpg(self, dir_basename=''):
        dirs = [dir_basename + 'head_met',
                dir_basename + 'head_unmet',
                dir_basename + 'tail_met',
                dir_basename + 'tail_unmet']
        coding_coords = self.__coding_coordinate()
        raw_pos = [self._get_raw_pos(i) for i in coding_coords]
        cpg_table = self.guppy_mbt[:,3]
        head_met_ratio, tail_met_ratio = self.get_coding_methylation_pros('cpg')
        if head_met_ratio < 0.015:
            self._truncate_fast5_file(dirs[1] + '/' + self.read_id + '.fast5', 
                    raw_pos[0], raw_pos[1], read_id=self.read_id[:-2]+'01')
        elif head_met_ratio > 0.06:
            self._truncate_fast5_file(dirs[0] + '/' + self.read_id + '.fast5', 
                                      raw_pos[0], raw_pos[1], read_id=self.read_id[:-2]+'01')
        if tail_met_ratio < 0.015:
            self._truncate_fast5_file(dirs[3] + '/' + self.read_id + '.fast5', 
                                      raw_pos[2], raw_pos[3], read_id=self.read_id[:-2]+'02')
        elif tail_met_ratio > 0.06:
            self._truncate_fast5_file(dirs[2] + '/' + self.read_id + '.fast5', 
                                      raw_pos[2], raw_pos[3], read_id=self.read_id[:-2]+'02')

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


def met_level_per_coding(reads, met_type, savename='temp.png'):
    met_pros = []
    for r in reads:
        head, tail = r.get_coding_methylation_pros(met_type)
        met_pros.append(head)
        met_pros.append(tail)
    weights = np.ones_like(met_pros)/float(len(met_pros))
    plt.hist(met_pros, bins=60, weights=weights)
    plt.savefig(savename, dpi=300)
    plt.close()


def plot_dam_positions(reads):
    pos_in_ref = []
    for r in reads:
        cpg_pos, dam_pos = r.get_methylated_bases(180)
        for pos in dam_pos:
            ref_pos = r.get_position_in_ref(pos, ref)
            if ref_pos:
                pos_in_ref.append(ref_pos)
    with open('dam_distribution.txt', 'w') as fw:
        fw.write(str(sorted(pos_in_ref)))
    plt.hist(pos_in_ref, bins=1000)
    plt.savefig('dam_distribution.png', dpi=500)
    return pos_in_ref


def basecalling(dirname, out_name):
    subprocess.run('~/Softwares/ont-guppy_4.2.2_linux64/ont-guppy/bin/guppy_basecaller -i ' + dirname + ' -s ' + out_name + ' -c dna_r9.4.1_450bps_modbases_dam-dcm-cpg_hac_prom.cfg --device cuda:0 --fast5_out', shell=True)


def resquiggle(bc_dirname, ref):
    subprocess.run('tombo resquiggle ' + bc_dirname + ' ' + ref + ' --processes 6 --num-most-common-errors 5', shell = True)


def tombo_plot(coordinates, offset, bc_dirname, bc_alt_dirname, pdf_filename):
    """
    coordinates, bc_dirname, bc_alt_dirname should be list
    """
    coord_str = ''
    for co in (np.array(coordinates) + offset) % 42999:
        coord_str += ' gi\|555853\|gb\|U13369.1\|HSU13369:' + str(co)
    bc_str = ''
    for dirname in bc_dirname:
        bc_str += dirname + ' '
    bc_alt_str = ''
    for dirname in bc_alt_dirname:
        bc_alt_str += dirname + ' '
    subprocess.run('tombo plot genome_locations --plot-standard-model --fast5-basedirs ' + bc_str + '--control-fast5-basedirs ' + bc_alt_str + '--genome-locations' + coord_str + ' --pdf-filename ' + pdf_filename, shell=True)


def separate_by_dam_and_plot_main(reads, coordinate, offset, ref):
    fast5_dir = 'dam_fast5s'
    bc_dir = fast5_dir + '_basecalled'
    control_dir = 'dam_control_fast5s'
    ctl_bc_dir = control_dir + '_basecalled'
    if os.path.exists(fast5_dir):
        shutil.rmtree(fast5_dir)
    if os.path.exists(control_dir):
        shutil.rmtree(control_dir)
    if os.path.exists(bc_dir):
        shutil.rmtree(bc_dir)
    if os.path.exists(ctl_bc_dir):
        shutil.rmtree(ctl_bc_dir)
    os.mkdir(fast5_dir)
    os.mkdir(control_dir)
    for r in reads:
        r.truncate_and_separate_by_dam(coordinate - 5, coordinate + 5, fast5_dir, control_dir, offset)
    basecalling(fast5_dir, bc_dir)
    basecalling(control_dir, ctl_bc_dir)
    shutil.rmtree(fast5_dir)
    shutil.rmtree(control_dir)
    resquiggle(bc_dir, ref)
    resquiggle(ctl_bc_dir, ref)
    tombo_plot(coordinate, offset, bc_dir, ctl_bc_dir, 
            pdf_filename='dam_' + str(coordinate) + '.pdf')


def separate_by_cpg_align_and_squiggle(reads, dir_basename=''):
    dirs = [dir_basename + 'head_met',
            dir_basename + 'head_unmet',
            dir_basename + 'tail_met',
            dir_basename + 'tail_unmet']
    called_dirs = [i + '_basecalled' for i in dirs]
    for directory in dirs:
        if os.path.exists(directory):
            shutil.rmtree(directory)
    for directory in dirs:
        os.mkdir(directory)
    for r in reads:
        r.truncate_and_separation_by_cpg(dir_basename)
    for in_dir, out_dir in zip(dirs, called_dirs):
        basecalling(in_dir, out_dir)
        shutil.rmtree(in_dir)
    for called_dir in called_dirs:
        resquiggle(called_dir, ref)


def plot_by_cpg(reads, coordinates, offset, ref, make_fast5s, dir_basename=''):
    dirs = [dir_basename + 'head_met',
            dir_basename + 'head_unmet',
            dir_basename + 'tail_met',
            dir_basename + 'tail_unmet']
    called_dirs = [i + '_basecalled' for i in dirs]
    if make_fast5s == 'yes':
        separate_by_cpg_align_and_squiggle(reads, dir_basename)
    bc_dirs = [i for i in called_dirs if '_met' in i]
    ctl_bc_dirs = [i for i in called_dirs if '_unmet' in i]
    coordinates = np.array(coordinates)
    tombo_plot(coordinates, bc_dirs, ctl_bc_dirs, 
            pdf_filename='1020cpg_' + str(coordinates[0]) + '-' + str(coordinates[-1]) + '.pdf')


def extract_rDNA_and_make_reads(rDNA_name_file, fast5_dir):
    if fast5_dir[-1] != '/':
        fast5_dir = fast5_dir + '/'
    rDNA_read_ids = []
    with open(rDNA_name_file) as f:
        for line in f:
            rDNA_read_ids.append(line.split()[0])
    fast5files = glob.glob(fast5_dir)
    rDNA_files = [fast5_dir+ i + '.fast5' for i in rDNA_read_ids]
    reads = []
    n = 0
    for rfile in rDNA_files:
        reads.append(Read(rfile, ref))
    return reads


def nanopore_reads_to_fig(in_dir, base_name, ref):
    if os.path.exists(base_name + '_plot'):
        shutil.rmtree(base_name + '_plot')
    os.mkdir(base_name + '_plot')
    subprocess.run('multi_to_single_fast5 -t 6 -i in_dir -s ' + base_name + '_single', shell=True)
    basecalling(base_name + '_single', base_name + '_bc_fast5s')
    shutil.rmtree(base_name + '_single')
    fast5s = glob.glob(base_name + '_bc_fast5s/workspace/*.fast5')
    ret_str = ''
    for f in fast5s:
        read = Read(fast5, ref)
        if read._is_this_healthy_rDNA():
            ret_str += '\n'.join(read.fastq) + '\n'
            read.plot_structure(ref, base_name + '_plot')
    with open(base_name + '_rDNAs.fastq', 'w') as fw:
        fw.write(ret_str)


if __name__ == '__main__':
    """
    CpG plotting requires a list or tuple of coordinates.
    dam plotting requires a scalar coordinate.
    """
    ref = 'rDNA_index/rDNA_for_cas9.fa'
    ref = 'rDNA_index/humRibosomal.fa'
    offset = get_ref_offset(ref)
    reads = extract_rDNA_and_make_reads('1026_rDNA_reads.txt', '1026_60_bc_fast5s/workspace/')
    met_level_per_coding(reads, 'cpg', '1026.png')
    quit()
    separate_by_cpg_align_and_squiggle(reads, 'rpa_')
    dirs1 = ['1015_head_unmet_basecalled/']
    dirs2 = ['1020_head_unmet_basecalled/']
    tombo_plot([420, 12000, 800], offset, dirs1, dirs2, 'comparison.pdf')
    quit()
    reads = extract_rDNA_and_make_reads('1015_dir_met_summary.txt', '1015_mod_fast5s/workspace/')
    separate_by_cpg_align_and_squiggle(reads, '1015_')
    intervals = range(11, 301, 20)
    plot_by_cpg(reads, intervals, offset, ref, make_fast5s='no', dir_basename='1020')
    result = re.finditer(r'GATC', standard_ref_seq)
    for match in result:
        gatc_pos = match.span()[0] + 2
        if 50 < gatc_pos < 13200:
            separate_by_dam_and_plot_main(reads, gatc_pos, offset, ref)
    quit()
