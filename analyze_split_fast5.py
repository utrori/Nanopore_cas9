import h5py
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


class Read(object):
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

    def __init__(self, fast5file, ref):
        self.fast5path = fast5file
        with h5py.File(fast5file) as f:
            read_name = list(f['Raw/Reads'].keys())[0]
            self.read_id = f['Raw/Reads/' + read_name].attrs['read_id'].decode()
            self.raw_duration = f['Raw/Reads/' + read_name].attrs['duration']
            self.raw = f['Raw/Reads/'+ read_name + '/Signal'][:]
            # it's important to check the location of ModBaseProbs in Analyses
            # in case it is the second time basecalling
            fastq = f['/Analyses/Basecall_1D_000/BaseCalled_template/Fastq']\
                    [()].decode().strip().split('\n')
            self.seq = fastq[1]
            self.quality = fastq[3]
            #fastq length! not raw signal length
            self.length = len(self.seq)
            self.guppy_mbt = \
                    f['/Analyses/Basecall_1D_000/BaseCalled_template/ModBaseProbs'][:]
            sam_summary = np.array(utilities.split_mapping_and_sam_analysis(300, 'temp', 
                    self.seq, self.quality, ref))
            self.ref_pos = sam_summary[:,2]
            self.direction = self._get_direction(sam_summary[:,1])

    def get_position_in_ref(self, pos, ref):
        with open('temp_files/temp_fastq.fastq', 'w') as fw:
            fw.write('@temp\n' + self.seq + '\n+\n' + 
                     self.quality)
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
                if current == pos:
                    ref_coordinate = ref_pos + pos - 1 + shift
                    break
            return ref_coordinate


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


def get_dam_pos(offset)
    result = re.finditer(r'GATC', standard_ref_seq)
    dam_pos = []
    for match in result:
        gatc_pos = match.span()[0] + 2
        if 50 < gatc_pos < 13200:
            dam_pos.append(gatc_pos + offset)
    return dam_pos


def get_and_copy_dam_file(coordinate, directory, ref, offset):
    fast5s = glob.glob(directory)
    pos_dir = 'pos_dir'
    neg_dir = 'neg_dir'


if __name__ == '__main__':
    ref = 'rDNA_index/rDNA_for_cas9.fa'
    offset = get_ref_offset(ref)
    r = Read('dam_control_fast5s_basecalled/workspace/0083d76f-a9eb-4662-ad07-093b0ce7db99.fast5')
