# -*- coding: utf-8 -*-

from search_rDNA_reads import make_temp_fastq
import subprocess
import os
import pandas as pd


def easy_flag(flag, base):
    """Return flag value based on the base.

    Args:
        flag (int): value of flag
        base (int): base
    Returns:
        the 1 or 0 value for the base value
    """
    return flag % (2 * base) // base


def analyze_split_reads(samfile, split_length):
    """Make transcribed region summary.

    Args:
        samfile (str): .sam filename
        split_length (int): split length
    Returns:
        TRs: a list of TR. A TR is a list of split mapped reads in the form of
        (position of split read, start coordinate in rDNA,
         end coordinate in rDNA, strand)
        When you have a continuous reads that are within transcribed region,
        they are grouped and named TR. Each TR should be each transcribed
        region of rDNA.
    """
    # this function receives a bwa output .sam file and output TRs.
    rDNA_len = 42999
    rDNA_TR_len = 13314
    with open(samfile) as sf:
        samdata = sf.readlines()[2:]
        """
        in_TR is a flag used to keep track of the where the split read is
        (is it inside TR or not). Not to confuse with flag (SAM header flag).
        TRs is a colletcion of read sets.
        Each read set is composed of reads that are mapped to each TR.
        """
        TRs = []
        TR = []
        prev_direction = ''
        TR_threshold = 10
        for n, item in enumerate(samdata):
            row = item.split()
            flag = int(row[1])
            mapped_coord = int(row[3])
            # reads that do not map to rDNA
            if easy_flag(flag, 4):
                continue
            # reads that are mapped to multiple loci are disregarded here
            if easy_flag(flag, 256):
                continue
            if easy_flag(flag, 16):
                direction = '-'
            else:
                direction = '+'
            if direction == '+':
                position = [mapped_coord,
                            (mapped_coord + split_length) % rDNA_len]
            else:
                position = [(mapped_coord - split_length) % rDNA_len,
                            mapped_coord]
            if 0 < position[0] < rDNA_TR_len and 0 < position[1] < rDNA_TR_len:
                # when the reads change its direction
                if prev_direction != '' and prev_direction != direction and \
                   len(TR) > TR_threshold:
                    TRs.append(TR)
                    TR = []
                # TR definition
                TR.append((n, position[0], position[1], direction))
            else:
                if len(TR) > TR_threshold:
                    TRs.append(TR)
                TR = []
            prev_direction = direction
    # end position to make the downstream analysis eaiser
    TRs.append([(n, 0, 0, '*')])
    return TRs


def analyze_all_fastq_file(fastq, split_length, length_cutoff):
    """Perform TR analysis for a fastq and output a summarized list.

    For each element in set_of_TRs, true TRs is element[1]!

    Args:
        fastq (str): fastq filename
        split_length (int): split length
        length_cutoff (int): reads shorter than this are omitted
    Returns:
        set_of_TRs: each TRs obtained from analyze_split_reads are paired with
        header as [header, TRs] and summarized as set_of_TRs
    """
    # this function read a fastq file and split each read and then map them to rDNA.
    # Each mapped read is analyzed by analyze_split_reads to make TRs.
    count = 0
    set_of_TRs = []
    with open(fastq) as f:
        for n, line in enumerate(f):
            if n % 4 == 0:
                header = line.strip()
            if n % 4 == 1:
                read = line.strip()
            if n % 4 == 3:
                quality = line.strip()
                if len(quality) < length_cutoff:
                    continue
                else:
                    count += 1
                    make_temp_fastq(split_length, header, read, quality)
                    FNULL = open(os.devnull, 'w')
                    subprocess.run('bwa mem -M -x ont2d -t 5 /home/yutaro/nanopore/clive/rDNA_index/humRibosomal.fa temp_files/temp_fastq.fastq > temp_files/temp_sam.sam', shell=True, stdout=FNULL, stderr=subprocess.STDOUT)
                    TRs = analyze_split_reads('temp_files/temp_sam.sam', split_length)
                    set_of_TRs.append([header, TRs])
    return set_of_TRs


def analyze_direction_distance(TRs, split_length):
    """Calculate direcion and distance between TRs.

    Args:
        TRs (list): TRs
        split_length (int): split length
    Returns:
        result: list of (direction, distance) between TRs
    """
    result = []
    for n in range(len(TRs) - 2):
        if TRs[n][0][3] == TRs[n+1][0][3]:
            direction = '='
        else:
            direction = '!'
        distance = split_length * (TRs[n+1][-1][0] - TRs[n][-1][0])
        result.append((direction, distance))
    return result


def visualize_TRs(TRs, split_length):
    """Output text visulaization of each TRs.

    Args:
        TRs (list): TRs
        split_length (int): split length
    Returns:
        output: text visualization of a read
    """
    size_of_bin = 1000
    prev_end = 0
    TR_lens = []
    dists = []
    directions = []
    prev_direction = ''
    output = ''
    if TRs == []:
        return output
    for TR in TRs[:-1]:
        # start and end are the numbers of split
        start = TR[0][0]
        direction = TR[0][3]
        directions.append(direction)
        end = TR[-1][0]
        TR_lens.append((end - start) * split_length)
        dists.append((start - prev_end) * split_length)
        prev_end = end
    for TR_len, dist, direction in zip(TR_lens, dists, directions):
        if direction == '+':
            output += '-' * int(dist / size_of_bin) + \
                      '￫' * int(TR_len / size_of_bin)
        else:
            output += '-' * int(dist / size_of_bin) + \
                      '￩' * int(TR_len / size_of_bin)
    dist = (TRs[-1][0][0] - prev_end) * split_length
    output += '-' * int(dist / size_of_bin)
    return output


def count_direction_dist(set_of_TRs, split_length):
    """Write text visualization of read for all TRs.

    Args:
        set_of_TRs (list): output of analyze_all_fastq_file
        split_length (int): split length
    Returns:
        None
    """
    summary = []
    for TRs in set_of_TRs:
        summary.extend(analyze_direction_distance(TRs[1], split_length))
    same = 0
    same_distances = []
    inverse = 0
    inverse_distances = []
    for item in summary:
        dist = item[1]
        if item[0] == '=':
            same += 1
            same_distances.append(dist)
        else:
            inverse += 1
            inverse_distances.append(dist)
    count = 0
    for dist in same_distances:
        if 35000 < dist < 50000:
            continue
        else:
            count += 1
    print(count)
    print((same, inverse))
    print(inverse_distances)


def analyze_fastq_by_header(fastq, header_id, split_length):
    """Find a read from a fastq file and analyze its structure.

    Args:
        fastq (str): filename
        header_id (str): header id
        split_length (int): split length
    Returns:
        TRs: TRs
    """
    with open(fastq) as f:
        for line in f:
            if header_id in line:
                header = line
                read = f.readline()
                f.readline()
                quality = f.readline()
                make_temp_fastq(split_length, header, read, quality)
                FNULL = open(os.devnull, 'w')
                subprocess.run('bwa mem -M -x ont2d -t 5 /home/yutaro/nanopore/'
                               'clive/rDNA_index/humRibosomal.fa '
                               'temp_files/temp_fastq.fastq > temp_files/'
                               'single_split_mapped.sam',
                               shell=True, stdout=FNULL,
                               stderr=subprocess.STDOUT)
                break
        TRs = analyze_split_reads('temp_files/single_split_mapped.sam', split_length)
        return TRs


if __name__ == '__main__':
    split_length = 200
    length_cutoff = 5000
    # this set_of_TRs is a set of [header, TRs]
    # TRs is a list of split reads that are mapped to the Trancribed Region. The number of the last split reads is added to the end of it. (see line 59)
    set_of_TRs = analyze_all_fastq_file('rDNA_containing_reads.fastq', split_length, length_cutoff)
    pd.to_pickle(set_of_TRs, 'set_of_TRs.pkl')
    #set_of_TRs = pd.read_pickle('set_of_TRs.pkl')
    summary = []
    for i in set_of_TRs:
        summary.append(analyze_direction_distance(i[1], split_length))
    for i in summary:
        if i:
            for item in i:
                if item[1] < 4000:
                    print(i)
    count_direction_dist(set_of_TRs, split_length)
    with open('reads_visualized.txt', 'w') as fw:
        for n, TRs in enumerate(set_of_TRs):
            fw.write(str(n) + '\t' + TRs[0].split()[0] + '\n' +  visualize_TRs(TRs[1], split_length) + '\n')
