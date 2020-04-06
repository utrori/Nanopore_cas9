import subprocess
import sys
import os
import itertools

FNULL = open(os.devnull, 'w')


def get_alingment_scores_in_sam(samfile, split_length):
    """Obtain alignment scores from split mapped samfile.

    The alignment scores are normalized by read length.

    Args:
        samfile (str): filename
    Returns:
        alignment_scores: a list of alingment scores of split reads
    """
    with open(samfile) as f:
        samdata = f.readlines()[2:]
        alignment_scores = []
        for line in samdata:
            row = line.split()
            """
            the second row of SAM file header is mapping status.
            256 means the mapping is shorter split mapping, not primary.
            the third row is contig on which the read is mapped.
            If the read is not mapped to any part of the index sequence,
            it shows *.
            """
            if row[2] != '*' and int(int(row[1]) % 512 / 256) != 1:
                read_length = len(row[9])
                """
                when the read_length is too short,
                it can produce very high score by chance
                """
                if read_length > 500:
                    alingment_score = int(line.split('AS:i:')[1].split()[0])
                    alignment_scores.append(alingment_score / split_length)
    return alignment_scores


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


def make_temp_fastq(split_length, header, read, quality,
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
    with open('temp_files/temp_fastq.fastq', 'w') as fw:
        for i in range(len(split_reads)):
            split_header = [header.split()[0], ' '.join(header.split()[1:])]
            fw.write(split_header[0] + '_' + str(i+1) + ' ' + split_header[1] +
                     '\n' + split_reads[i] + '\n+\n' +
                     split_qualities[i] + '\n')


def check_fastq_file(fastq_file, split_length):
    """Check reads in a fastq file and collect rDNA containing reads.

    Args:
        fastq_file (str): fastq filename
        split_length (int): split length
    Returns:
        None
    """
    with open(fastq_file) as f:
        for data in itertools.zip_longest(*[iter(f)]*4):
            header = data[0].strip()
            read = data[1].strip()
            quality = data[3].strip()
            make_temp_fastq(split_length, header, read, quality)
            """
            You can avoid printing the output of subprocess by directing
            the stdout to devnull.
            """
            subprocess.run('bwa mem -M -x ont2d -t 5 '
                           '/home/yutaro/nanopore/clive/'
                           'rDNA_index/humRibosomal.fa temp_files/temp_fastq.fastq > '
                           'temp_sam.sam', shell=True, stdout=FNULL,
                           stderr=subprocess.STDOUT)
            alignment_scores = get_alingment_scores_in_sam('temp_sam.sam', split_length)
            if any(AS > 0.5 for AS in alignment_scores):
                print(alignment_scores)
                print(len(read))
                with open('rDNA_containing_reads.fastq', 'a') as fw:
                    fw.write(header + '\n' + read + '\n+\n' + quality + '\n')


if __name__ == '__main__':
    """
    This process split each read into the size of split_length
    and map them to rDNA consensus sequence.
    """
    split_length = 4000
    """
    input_fastq = ('/home/yutaro/data/cliveome/fastq_runid_'
                   '0ee57c1a265c2f494821757929f1af60fe060a43_0'
                   '.fastq')
    """
    input_fastq = sys.argv[1]
    check_fastq_file(input_fastq, split_length)
