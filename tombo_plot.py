import subprocess
import numpy as np

r = range(10000, 23000, 20)
#r = [11803]
positions = ''.join([' gi\|555853\|gb\|U13369.1\|HSU13369:' + str(s) for s in r])
string = 'tombo plot genome_locations --genome-locations' + positions
string += ' --fast5-basedirs rpa_head_unmet_basecalled/ rpa_tail_unmet_basecalled/ --control-fast5-basedirs rpa_met_basecalled/ rpa_head_met_basecalled/ --plot-standard-model --pdf-filename rpa_2.pdf'
subprocess.run(string, shell=True)
quit()
string = 'tombo plot genome_locations --genome-locations' + positions
string += ' --fast5-basedirs rpa_head_unmet_basecalled/ rpa_tail_unmet_basecalled/ --plot-standard-model --plot-alternate-model 6mA --pdf-filename rpa_dam.pdf'
