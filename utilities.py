import subprocess
import os
import numpy as np
import matplotlib.collections as mc
import matplotlib.cm as cm
from matplotlib import pyplot as plt
import sys


FNULL = open(os.devnull, 'w')


def easy_flag(flag, base):
    """Return flag value based on the base.

    Args:
        flag (int): value of flag
        base (int): base
    Returns:
        the 1 or 0 value for the base value
    """
    return flag % (2 * base) // base


def circular_slice(sequence, index1, index2):
    """Return slice of a sequence even when the indexes are out of range.

    Note: this function does not work for indexes that are bigger than twice
    the size of sequence length.
    Index1 should be smaller than index2.

    Args:
        sequence (str): any sequence
        index1 (int): index for slicing
        index2 (int): index for slicing
    Returns:
        sliced_seq: sliced sequence
    """
    length = len(sequence)
    if 0 <= index1 <= length and 0 <= index2 <= length:
        return(sequence[index1:index2])
    if index1 < 0:
        return sequence[index1:] + circular_slice(sequence, 0, index2)
    if index2 > length:
        return sequence[index1:index2] + sequence[:index2 % length]


def split_sequence(sequence, split_length):
    """Split a fastq file based on the split length.

    Args:
        sequence (str): any string
        split_length (int): the string is split by the given length
    Returns:
        split_seq: a list of split strings
    """
    split_seq = []
    for i in range(int(len(sequence) / split_length + 1)):
        split_seq.append(sequence[i*split_length:(i+1)*split_length])
    return split_seq


def make_temp_fastq(split_length, header, read, quality, identifier
                    tempfilename='temp_files/temp_fastq.fastq'):
    """Make temporary fastq by splitting the read by split length.

    Args:
        split_length (int): split length
        header (str): header
        read (str): read sequence
        quality (str): quality track of fastq
        tempfilename (str): temp file name
    Returns:
        None
    """
    split_reads = split_sequence(read, split_length)
    split_qualities = split_sequence(quality, split_length)
    if header[0] != '@':
        header = '@' + header
    with open('temp_files_' + identifier + '/temp_fastq.fastq', 'w') as fw:
        for i in range(len(split_reads)):
            split_header = [header.split()[0], ' '.join(header.split()[1:])]
            fw.write(split_header[0] + '_' + str(i+1) + ' ' + split_header[1] +
                     '\n' + split_reads[i] + '\n+\n' +
                     split_qualities[i] + '\n')


def bwa_mapping(ref_path, in_fastq, out_sam, multi=False):
    if multi == True:
        subprocess.run('bwa mem -Ma -x ont2d -t 6 {} {} > {}'.format(ref_path, in_fastq, out_sam),
                   shell=True, stdout=FNULL, stderr=subprocess.STDOUT)
    else:
        subprocess.run('bwa mem -M -x ont2d -t 6 {} {} > {}'.format(ref_path, in_fastq, out_sam),
                   shell=True, stdout=FNULL, stderr=subprocess.STDOUT)


def split_mapping_and_sam_analysis(split_length, header, read, quality, ref):
    identifier = str(random.random())
    make_temp_fastq(split_length, header, read, quality, identifier)
    bwa_mapping(ref, 'temp_files/temp_fastq_' + identifier + '.fastq', 'temp_files/single_split_mapped_' + identifier + '.sam')
    sam_info = []
    with open('temp_files/single_split_mapped_' + identifier + '.sam') as f:
        samdata = f.readlines()[2:]
        split_info = []
        # split_info is like [flag, direction, position, CIGAR, score]
        for line in samdata:
            row = line.split()
            flag = int(row[1])
            split_info.append(flag)
            if easy_flag(flag, 4):
                split_info.append(0)
            elif easy_flag(flag, 16):
                split_info.append('-')
            else:
                split_info.append('+')
            split_info.append(int(row[3]))
            split_info.append(row[5])
            if flag != 4:
                split_info.append(int(line.split('\tAS:i:')[1].split()[0]))
            else:
                split_info.append(0)
            sam_info.append(split_info)
            split_info = []
    os.remove('temp_files/temp_fastq_' + identifier + '.fastq')
    os.remove('temp_files/single_split_mapped_' + identifier + '.sam')
    return np.array(sam_info)


def make_shifted_ref(ref_filename, new_ref_filename, shift):
    with open(new_ref_filename, 'w') as fw:
        with open(ref_filename) as f:
            header = f.readline()
            fw.write(header)
            seq = ''
            for line in f:
                seq += line.strip()
            new_seq = circular_slice(seq, shift, len(seq)+shift)
            print(len(new_seq))
            for item in [new_seq[i:i+70] for i in range(0, len(new_seq), 70)]:
                fw.write(item + '\n')


def plot_read_structure(header, split_length, offset=0, savename=None, title=None):
    """Plot read based on split mapped sam file.

    Args:
        header (str): header name
        split_length (int): split_length of the read when it was mapped to rDNA
        offset (int): offset of rDNA reference file
        savename (str): if this is specifed, the plot is saved. otherwise,
        it shows the plot.

    Returns:
        None
    """
    samfile = 'temp_files/single_split_mapped.sam'

    plot = []
    rD_size = 42999
    TR_size = 13314
    with open(samfile) as f:
        samdata = f.readlines()[2:]
        for line in samdata:
            row = line.split()
            flag = int(row[1])
            # temporarily disregard multiply mapped reads
            if easy_flag(flag, 4) == 1:
                plot.append((-10000, '*'))
            else:
                if easy_flag(flag, 16) != 1:
                    plot.append(((int(row[3]) + offset) % 42999, '+'))
                else:
                    plot.append(((int(row[3]) + offset) % 42999, '-'))
    read_num = len(plot)

    x = np.linspace(0, split_length * (read_num - 1), read_num)
    r_coord = [i[0] for i in plot]
    vertical_lines = []
    for n in range(read_num):
        if plot[n][1] == '+':
            line_element = [(x[n], r_coord[n]),
                            (x[n] + split_length, r_coord[n] + split_length)]
        elif plot[n][1] == '-':
            line_element = [(x[n], r_coord[n]),
                            (x[n] + split_length, r_coord[n] - split_length)]
        else:
            line_element = [(x[n], r_coord[n]),
                            (x[n] + split_length, r_coord[n])]
        vertical_lines.append(line_element)
    vertical_lines.append([[0, 0], [split_length * read_num, 0]])
    vertical_lines.append([[0, 0 + TR_size],
                          [split_length * read_num, 0 + TR_size]])
    vertical_lines.append([[0, 0 + rD_size],
                          [split_length * read_num, 0 + rD_size]])

    lw = [1 for i in range(read_num)]
    lw.extend([0.5 for i in range(3)])
    ls = ['solid' for i in range(read_num)]
    ls.extend(['dashed' for i in range(3)])
    cl = []
    for coord in r_coord:
        if coord == -10000:
            cl.append('black')
        else:
            cl.append(cm.hsv((coord - offset) / rD_size))
    cl.extend(['black' for i in range(3)])
    lc = mc.LineCollection(vertical_lines, linewidths=lw,
                           linestyles=ls, colors=cl)

    if savename:
        fig = plt.figure()
        plt.subplots_adjust(left=0.2)
        ax = fig.add_subplot()
        ax.add_collection(lc)
        ax.autoscale()
        ax.set_yticks((-10000, 0, 10000, 20000, 30000, 40000))
        ax.set_yticklabels(('unmapped', 0, 10000, 20000, 30000, 40000))
        ax.set_ylim((-12000, 45000))
        if title:
            ax.set_title(title)
        if not savename:
            plt.show()
        else:
            plt.savefig(savename, dpi=300)
        plt.close()
    else:
        return lc


if __name__ == '__main__':
    make_shifted_ref('rDNA_index/humRibosomal.fa', 'rDNA_index/rDNA_for_cas9.fa', -10000)
