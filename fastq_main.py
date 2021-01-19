import itertools 
import utilities
import numpy as np


ref = 'rDNA_index/humRibosomal.fa'
with open('wstemp.fastq') as f:
    for items in itertools.zip_longest(*[iter(f)]*4):
        if len(items[1]) > 10000:
            sam_temp = utilities.split_mapping_and_sam_analysis(300, items[0].split()[0][1:], items[1], items[3], ref)
            sam_info = np.array(sam_temp)
            print(sam_info)
