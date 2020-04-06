from find_reads_at_the_ends import circular_slice

rDNA = ''
with open('/home/yutaro/nanopore/clive/rDNA_index/humRibosomal.fa') as f:
    f.readline()
    for line in f:
        rDNA += line.strip()

# it is better to have some margin at the both ends
margin = 2000
offset = 9000
sequence = circular_slice(rDNA, offset, 42999 + offset + margin)

with open('reference/rDNA_for_cas9.fasta', 'w') as fw:
    fw.write('>rDNA\n' + sequence)
