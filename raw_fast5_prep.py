import search_rDNA_reads
import os
import shutil
import subprocess
import glob

filename = '/var/lib/minknow/data/200825_60_cas9/60/20200825_0808_MN32877_ADX215_125bb26f/fast5_pass/'
name_pref = '0825_60'

if os.path.exists(name_pref):
    shutil.rmtree(name_pref)
if os.path.exists(name_pref + '_single'):
    shutil.rmtree(name_pref + '_single')
if os.path.exists(name_pref + '_bc_fast5s'):
    shutil.rmtree(name_pref + '_bc_fast5s')
shutil.copytree(filename, name_pref, symlinks=False, ignore=None)
subprocess.run('multi_to_single_fast5 -t 5 -i ' + name_pref + ' -s ' + name_pref + '_single', shell=True)
subprocess.run('~/Softwares/ont-guppy_4.2.2_linux64/ont-guppy/bin/guppy_basecaller -i ' + name_pref + '_single -s ' + name_pref + '_bc_fast5s -c dna_r9.4.1_450bps_modbases_dam-dcm-cpg_hac_prom.cfg --device cuda:0 --fast5_out --recursive', shell=True)
fast5s = glob.glob(name_pref + '_bc_fast5s/workspace/**/*.fast5')
for f in fast5s:
    shutil.move(f, name_pref + '_bc_fast5s/workspace/')
search_rDNA_reads.search_rDNA_reads(name_pref + '_bc_fast5s', name_pref + '_rDNA_reads.txt')
shutil.rmtree(name_pref)
shutil.rmtree(name_pref + '_single')
