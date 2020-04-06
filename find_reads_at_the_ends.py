from search_rDNA_reads import make_temp_fastq
from rDNA_structure_of_long_reads import easy_flag
from visualize_nanopore_read import plot_read_structure
from Bio.Seq import Seq
import os
import subprocess
import itertools
import pandas as pd
import matplotlib.pyplot as plt


FNULL = open(os.devnull, 'w')


def read_continuity(coordinates, split_length):
    """Judge if the list of coordinates are cotinuous.
    
    Args:
        coordinates (list(int, int,...)): list of coordinates
        split_length (int): split_length
    Returns:
        1 or 0
    """
    for i in range(len(coordinates) - 1):
        if split_length - 100 < abs(coordinates[i+1] - coordinates[i]) \
           < split_length + 100:
            continue
        else:
            return 0
    return 1


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


def find_end_reads(samfile, split_length, rDNA_coordinate=0):
    """Check if the read is at the end of rDNA based on split mapped .sam.

    This function returns 0 if the read is not at the end. If it is at the rDNA
    end, this returns rough information of the rDNA end.

    Args:
        samfile (str): path of split mapped sam
        split_length (int): length of split when the read is mapped with bwa
        rDNA_coordinate (int): if 1, also returns rDNA coorinate

    Returns:
        0: if the read is not at the end
        split_length * n: provisional cooridate of rDNA end
        direction: the direciton of rDNA, '+' or '-'
        side_of_non_rDNA: 'left' or 'right'
        rDNA_coordiante: rDNA coordinate
    """
    with open(samfile) as f:
        samdata = f.readlines()[2:]
    mapping_state = []
    threshold = 0.1
    for item in samdata:
        row = item.split()
        flag = int(row[1])
        if easy_flag(flag, 4) != 1:
            if easy_flag(flag, 16) != 1:
                mapping_state.append(1)
            else:
                mapping_state.append(-1)
            # 1 for + direction and -1 for - direction
        else:
            mapping_state.append(0)

    offset = 70
    # you have to ignore the end of reads
    # because it can cause very high non-mapping rate by chance
    candidates = []
    for n in range(1, len(mapping_state)-1):
        left = [abs(i) for i in mapping_state[:n]]
        right = [abs(i) for i in mapping_state[n:]]
        left_rate = sum(left) / len(left)
        right_rate = sum(right) / len(right)
        if left_rate > (1 - threshold) and right_rate < threshold:
            candidates.append((n, left_rate, right_rate))
        elif left_rate < threshold and right_rate > (1 - threshold):
            candidates.append((n, left_rate, right_rate))

    if not candidates:
        return 0
    else:
        provisional_best_diff = 0
        provisional_best_pos = 0
        for n, lr, rr in candidates:
            diff = abs(lr - rr)
            if diff > provisional_best_diff:
                provisional_best_diff = diff
                provisional_best_pos = n
                provisional_lr = lr
                provisional_rr = rr
        lr = provisional_lr
        rr = provisional_rr
        n = provisional_best_pos

        if n <= offset or (len(mapping_state) - n) <= offset:
            return 0

        if lr > rr:
            side_of_non_rDNA = 'right'
            near_boundary_mss = [i for i in mapping_state[n-20:n]]
            rDNA_coordinate = int(samdata[n-1].split()[3])
        else:
            side_of_non_rDNA = 'left'
            near_boundary_mss = [i for i in mapping_state[n:n+20]]
            rDNA_coordinate = int(samdata[n].split()[3])
            return 0
        # average direction of 20 reads at the boudary region
        if sum(near_boundary_mss) > 0:
            direction = '+'
        else:
            direction = '-'
        if rDNA_coordinate:
            return (split_length * n, direction, side_of_non_rDNA,
                    rDNA_coordinate)
        else:
            return (split_length * n, direction, side_of_non_rDNA)


def make_fastq_for_boundary(fastq, boundary, half_length):
    """Make fastq for the determination of exact coordinate of rDNA end.

    Args:
        fastq (list(str, str, str)): (header, read, quality)
        boundary (int): Rough coordinate of boundary inferred by find_end_reads
        half_length (int): half length of temp fastq

    Returns:
        Nothing. This script generates temp_boundary.fastq file.
    """
    half_length = 1000
    output = fastq[0] + '\n'
    output += fastq[1][boundary - half_length:boundary + half_length] + '\n'
    output += '+' + '\n'
    output += fastq[2][boundary - half_length:boundary + half_length]
    with open('temp_bound.fastq', 'w') as fw:
        fw.write(output)


def find_true_boundary(header, read, quality, temp_boundary):
    """Map fastq for boundary and determine the boundary coordinate.

    Args:
        header (str): header
        read (str): read sequence
        quality (str): read quality
        temp_boundary (list(int, str, str)): (rough boundary coordinate,
                                              direction, side_of_non_rDNA)

    Returns:
        true_boundary: boundary coordinate. '' if the result is ambiguous.
        rDNA_boundary: boundary in terms of rDNA coordinate
    """
    half_length = 1000
    make_fastq_for_boundary((header, read, quality),
                            temp_boundary[0], half_length)
    subprocess.run('bwa mem -M -x ont2d -t 5 '
                   '/home/yutaro/nanopore/clive/rDNA_index/humRibosomal.fa '
                   'temp_files/temp_bound.fastq > temp_files/temp_bound.sam',
                   shell=True, stdout=FNULL, stderr=subprocess.STDOUT)
    with open('temp_files/temp_bound.sam') as samf:
        samdata = samf.readlines()[2:]
    row = samdata[0].split()
    rDNA_coordinate = int(row[3])
    direction = temp_boundary[1]
    side_of_non_rDNA = temp_boundary[2]
    cigar = row[5]
    # only when the CIGAR sequence is simple
    true_boundary = ''
    rDNA_boundary = ''
    if side_of_non_rDNA == 'right':
        temp_dist = cigar.split('M')[-1]
        if 'S' in temp_dist:
            dist_to_boundary = int(temp_dist[:-1])
            # non rDNA reads are expressed as read[true_boundary:]
            true_boundary = temp_boundary[0] + half_length - dist_to_boundary
            if direction == '+':
                rDNA_boundary = (rDNA_coordinate + 2 * half_length
                                 - dist_to_boundary) % 42999
            else:
                rDNA_boundary = (rDNA_coordinate - 2 * half_length
                                 + dist_to_boundary) % 42999
    else:
        temp_dist = cigar.split('S')[0]
        if temp_dist.isdigit():
            dist_to_boundary = int(temp_dist)
            # non rDNA reads are expressed as read[:true_boundary]
            true_boundary = temp_boundary[0] - half_length + dist_to_boundary
            if direction == '+':
                rDNA_boundary = (rDNA_coordinate + dist_to_boundary) % 42999
            else:
                rDNA_boundary = (rDNA_coordinate - dist_to_boundary) % 42999
    return (true_boundary, rDNA_boundary)


def find_true_boundary2(header, read, quality, temp_boundary):
    """Map fastq for boundary and determine the boundary coordinate.

    In version 2, rough rDNA position estimate by the find_end_reads
    is also used.

    Args:
        header (str): header
        read (str): read sequecne
        quality (str): read quality
        temp_boundary (list): (rough boundary, direction, side_of_non_rDNA,
                               rough rDNA boundary)
    Returns:
        true_boundary: boundary coordinate. '' if the result is ambiguous.
        rDNA_boundary: boundary in terms of rDNA coordinate
    """
    rDNA_seq = ''
    with open('rDNA_index/humRibosomal.fa') as f:
        fadata = f.readlines()
    for line in fadata[1:]:
        rDNA_seq += line.strip()

    index_half_len = 2000
    # r_bound is rough boundary in terms of read
    r_bound = temp_boundary[0]
    temp_index = read[r_bound-index_half_len:r_bound+index_half_len]
    with open('temp_index/temp_index.fasta', 'w') as fw:
        fw.write('>temp_index\n' + temp_index)
    subprocess.run('bwa index temp_index/temp_index.fasta', shell=True,
                   stdout=FNULL, stderr=subprocess.STDOUT)

    rD_half = 2000
    # rD_bound is rought boundary in terms of rDNA coordiante
    rD_bound = temp_boundary[3]
    parital_rDNA = circular_slice(rDNA_seq, rD_bound-rD_half, rD_bound+rD_half)
    with open('temp_files/temp_rDNA.fastq', 'w') as fw:
        fw.write('>temp\n' + parital_rDNA + '\n+\n' + 'J' * 2 * rD_half)

    subprocess.run('bwa mem -M -x ont2d -t 5 '

                   '/home/yutaro/nanopore/clive/temp_index/temp_index.fasta '
                   'temp_files/temp_rDNA.fastq > temp_files/temp_bound.sam',
                   shell=True, stdout=FNULL, stderr=subprocess.STDOUT)
    with open('temp_files/temp_bound.sam') as samf:
        samdata = samf.readlines()[2:]
    row = samdata[0].split()
    coordinate = int(row[3])
    direction = temp_boundary[1]
    side_of_non_rDNA = temp_boundary[2]
    cigar = row[5]

    if side_of_non_rDNA == 'right':
        temp_dist = cigar.split('M')[-1][:-1]
    else:
        temp_dist = cigar.split('S')[0]
    if temp_dist.isdigit():
        temp_dist = int(temp_dist)
        if temp_dist < 400:
            return 0
    else:
        return 0


def find_boundaries_from_fastq(fastq, split_length):
    """Find boundary containing reads from a .fastq file.

    Args:
        fastq (str): filename
        split_length (split_length): split_length
    Returns:
        boundaries: list of (boundary coordinate, boundary rDNA coordinate,
                             direction, side of non rDNA, read, header)
    """
    with open(fastq) as f:
        boundaries = []
        for n, each_fastq in enumerate(itertools.zip_longest(*[iter(f)]*4)):
            header = each_fastq[0].strip()
            read = each_fastq[1].strip()
            if len(read) < 40000:
                continue
            quality = each_fastq[3].strip()
            make_temp_fastq(split_length, header, read, quality)
            subprocess.run('bwa mem -M -x ont2d -t 5 '
                           '/home/yutaro/nanopore/clive/rDNA_index/'
                           'humRibosomal.fa temp_files/temp_fastq.fastq > '
                           'temp_files/temp_sam.sam',
                           shell=True, stdout=FNULL, stderr=subprocess.STDOUT)
            # rDNA_coordinate=1 when using find_true_boundary2
            temp_boundary = find_end_reads('temp_files/temp_sam.sam',
                                           split_length, rDNA_coordinate=1)
            if temp_boundary:
                header = header.split()[0]
                boundary = int(temp_boundary[0])
                direction = temp_boundary[1]
                side = temp_boundary[2]
                if direction == '+':
                    if side == 'right':
                        bound_seq = read[boundary:boundary+10000]
                        with open('boundary_seq1.fa', 'a') as fw:
                            fw.write('>' + header + '\n')
                            fw.write(bound_seq + '\n\n')
                    else:
                        with open('boundary_seq2.fa', 'a') as fw:
                            fw.write('>' + header + '\n')
                            fw.write(bound_seq + '\n\n')
                else:
                    bound_seq = read[boundary-10000:boundary]
                    revcom = str(Seq(bound_seq).reverse_complement())
                    if side == 'left':
                        with open('boundary_seq1.fa', 'a') as fw:
                            fw.write('>' + header + '\n')
                            fw.write(revcom + '\n\n')
                    else:
                        with open('boundary_seq2.fa', 'a') as fw:
                            fw.write('>' + header + '\n')
                            fw.write(revcom + '\n\n')

                plot_read_structure(header, split_length, savename='end_reads/' +
                                    header + '.png', title=str(temp_boundary[0]))
                continue
                true_boundary = find_true_boundary2(header, read,
                                                    quality, temp_boundary)
                if true_boundary:
                    boundaries.append((true_boundary[0], true_boundary[1],
                                       temp_boundary[1], temp_boundary[2],
                                       read, header.split()[0]))
        return boundaries


if __name__ == '__main__':
    fastq = 'rDNA_reads.fastq'
    split_length = 200

    boundaries = find_boundaries_from_fastq(fastq, split_length)
    pd.to_pickle(boundaries, 'pkls/rDNA_boundaries.pkl')
    #boundaries = pd.read_pickle('pkls/rDNA_boundaries.pkl')
    quit()
    rDNA_coords1 = []
    rDNA_coords2 = []
    for item in boundaries:
        if item[2] == '+':
            if item[3] == 'left':
                rDNA_coords1.append(item[1])
            else:
                rDNA_coords2.append(item[1])
        else:
            if item[3] == 'left':
                rDNA_coords2.append(item[1])
            else:
                rDNA_coords1.append(item[1])

    plt.subplot(2, 1, 1)
    plt.hist(rDNA_coords1, bins=50, range=(0, 45000))
    plt.subplot(2, 1, 2)
    plt.hist(rDNA_coords2, bins=50, range=(0, 45000))
    plt.show()
