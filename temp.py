import itertools

mega = {}
with open('megalodon_results/basecalls.fasta') as f:
    for item in itertools.zip_longest(*[iter(f)]*2):
        read_id = item[0][1:].strip()
        read = item[1].strip()
        mega[read_id] = read

guppy = {}
with open('FAL81148_pass_885db4b7_0.fastq') as f:
    for item in itertools.zip_longest(*[iter(f)]*4):
        read_id = item[0][1:].split()[0]
        read = item[1].strip()
        guppy[read_id] = read

rDNAs = []
with open('reads_visualized.txt') as f:
    for item in itertools.zip_longest(*[iter(f)]*2):
        rDNAs.append(item[0].strip().split()[1][1:])

for read_id in rDNAs:
    print(len(mega[read_id]))
    print(len(guppy[read_id]))
