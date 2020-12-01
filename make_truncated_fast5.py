import h5py
import numpy as np
import glob
import itertools
import matplotlib.pyplot as plt
import utilities
from hmmlearn import hmm


def rDNA_read_ids2seq(filename):
    ids2seq = {}
    with open(filename) as f:
        for items in itertools.zip_longest(*[iter(f)]*4):
            ids2seq[items[0].split()[0][1:]] = items[1].strip()
    return ids2seq


def get_direction(seq):
    summary = np.array((utilities.split_mapping_and_sam_analysis(200, 'test', seq, 'J' * len(seq), '/home/yutaro/nanopore/clive/rDNA_index/humRibosomal.fa')))
    directions = summary[:, 1]
    minus = 0.1
    plus = 0.1
    for i in directions:
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


def get_id2fast5(ids, fast5s):
    id2fast5 = {}
    for f5 in fast5s:
        fid = f5.split('/')[-1].split('.')[0]
        if fid in ids:
            id2fast5[fid] = f5
    return id2fast5


def coding_coordinate(direction, length):
    region1 = 10500
    region2 = 3500
    if direction == '+':
        a_s = 0
        a_e = region2
        b_s = length - region1
        b_e = length - 1
    elif direction == '-':
        a_s = 0
        a_e = region1
        b_s = length - region2
        b_e = length - 1
    return (a_s, a_e, b_s, b_e)


def truncate_fast5_file(filename, outname, start, end):
    length = end - start
    with h5py.File(filename, 'r') as f:
        with h5py.File(outname,'w') as fw:
            f.copy('Raw', fw)
            f.copy('UniqueGlobalKey', fw)
            read_name = list(f['Raw/Reads'].keys())[0]
            fw['/'].attrs.create('file_version', 2, dtype='double')
            fw['Raw/Reads/' + read_name].attrs['duration'] = length
            signal = f['Raw/Reads/'+ read_name + '/Signal']
            del fw['Raw/Reads/'+ read_name + '/Signal']
            temp = fw['Raw/Reads/' + read_name].create_dataset('Signal', (length,), dtype='<i2', maxshape=(None,), chunks=(length,), compression=1)
            temp[:] = signal[start:end]
            fw.flush()


def get_raw_pos(filename, coding_co, rlen):
    with h5py.File(filename, 'r') as f:
        read_name = list(f['Raw/Reads'].keys())[0]
        signal = np.array((f['Raw/Reads/' + read_name + '/Signal']))
        slen = len(signal)
        a_s = 0
        a_e = int(slen/rlen * coding_co[1])
        b_s = int(slen/rlen * coding_co[2])
        b_e = slen
    return (a_s, a_e, b_s, b_e)


def get_cpg_met(rid, filename):
    with open(filename) as f:
        for line in f:
            row = line.strip().split()
            if row[0] == rid:
                return (float(row[2]), float(row[3]))


def get_direction_from_file(rid, filename):
    with open(filename) as f:
        for line in f:
            row = line.strip().split()
            if row[0] == rid:
                return row[1]


if __name__ == '__main__':
    id2seq = rDNA_read_ids2seq('1015_rDNA_reads.fastq')
    fast5s = glob.glob('1015_single_fast5s/*/*.fast5')
    id2fast5 = get_id2fast5(list(id2seq.keys()), fast5s)
    for rid, seq in id2seq.items():
        filename = id2fast5[rid]
        length = len(seq)
        direction = get_direction_from_file(rid, '1015_dir_met_summary.txt')
        if not direction or length < 39000 or length > 46000:
            continue
        coding_coords = coding_coordinate(direction, length)
        raw_coords = get_raw_pos(filename, coding_coords, length)
        cpg_met = get_cpg_met(rid, '1015_dir_met_summary.txt')
        print(cpg_met)
        if cpg_met[0] < 15:
            truncate_fast5_file(filename, '1015_truncated_met1/' + rid + '_1.fast5', raw_coords[0], raw_coords[1])
        if cpg_met[1] < 15:
            truncate_fast5_file(filename, '1015_truncated_unmet1/' + rid + '_1.fast5', raw_coords[0], raw_coords[1])
        if cpg_met[0] > 35:
            truncate_fast5_file(filename, '1015_truncated_met2/' + rid + '_2.fast5', raw_coords[2], raw_coords[3])
        if cpg_met[1] > 35:
            truncate_fast5_file(filename, '1015_truncated_unmet2/' + rid + '_2.fast5', raw_coords[2], raw_coords[3])
        #truncate_fast5_file(filename, '1015_truncated1/' + rid + '_1.fast5', raw_coords[0], raw_coords[1])
        #truncate_fast5_file(filename, '1015_truncated2/' + rid + '_2.fast5', raw_coords[2], raw_coords[3])
