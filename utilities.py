import subprocess
import os
from search_rDNA_reads import make_temp_fastq
from rDNA_structure_of_long_reads import easy_flag

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


def split_mapping_and_sam_analysis(split_length, header, read, quality, ref):
    make_temp_fastq(split_length, header, read, quality)
    FNULL = open(os.devnull, 'w')
    subprocess.run('bwa mem -M -x ont2d -t 5 ' + ref + 
                   ' temp_files/temp_fastq.fastq > temp_files/'
                   'single_split_mapped.sam',
                   shell=True, stdout=FNULL,
                   stderr=subprocess.STDOUT)
    sam_info = []
    with open('temp_files/single_split_mapped.sam') as f:
        samdata = f.readlines()[2:]
        split_info = []
        # split_info is like [mapping_status, direction, position, CIGAR, score]
        for line in samdata:
            row = line.split()
            flag = int(row[1])
            if easy_flag(flag, 4):
                split_info.append(0)
            else:
                split_info.append(1)
            if easy_flag(flag, 16):
                split_info.append('-')
            else:
                split_info.append('+')
            split_info.append(int(row[3]))
            split_info.append(row[5])
            split_info.append(int(line.split('AS:i:')[1].split()[0]))
            sam_info.append(split_info)
            split_info = []
    return sam_info
