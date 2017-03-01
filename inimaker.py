from __future__ import print_function,division
import numpy as np
from os.path import join
import os
import argparse

def main(datadir,tag,param_sampler,*samples):
    cosmosis_root = '/share/splinter/hj/PhD/cosmosis'
    if samples[0][0]=='all':
	print('GENERATING .ini files for all high/low red/blue samples...!')
        samples = [['hzr','hzb','lzr','lzb']]
    for sample in samples[0]:
        if 'h' in sample:
            z = 'high'
        elif 'l' in sample:
            z = 'low'
        if 'r' in sample:
            col = 'Red'
        elif 'b' in sample:
            col = 'Blue'
        if param_sampler=='grid':
            iniparams = [
            '[runtime]\n',
            'sampler = grid\n',
            '\n',
            '[DEFAULT]\n',
            'N_BIN=50\n',
            'N_Z=50\n',
            'MY_ROOT=/share/splinter/hj/PhD\n',
            'DATA_PATH=%%(MY_ROOT)s/%s\n'%datadir,
            'TAG=%s\n'%tag,
            'Z_SAMPLE=%s\n'%z,
            'COLOUR=%s\n'%col,
            'SHAPEZ=%%(DATA_PATH)s/%%(Z_SAMPLE)sZ_%%(COLOUR)s_galZs.txt\n',
            'DENSZ=%%(DATA_PATH)s/%%(Z_SAMPLE)sZ_galZs.txt\n',
            'SAVE_NZ=%%(DATA_PATH)s/nofz\n',
            '\n',
            '[grid]\n',
            'fatal_errors=T\n',
            'nsample_dimension=10\n',

            '[output]\n',
            'format=text\n',
            'filename=IAs_%%(TAG)s_%%(Z_SAMPLE)sZ_%%(COLOUR)s\n'
            ]
        elif param_sampler=='test':
            iniparams = [
            '[runtime]\n',
            'sampler = test\n',
            '\n',
            '[DEFAULT]\n',
            'N_BIN=50\n',
            'N_Z=50\n',
            'MY_ROOT=/share/splinter/hj/PhD\n',
            'DATA_PATH=%%(MY_ROOT)s/%s\n'%datadir,
            'TAG=%s\n'%tag,
            'Z_SAMPLE=%s\n'%z,
            'COLOUR=%s\n'%col,
            'SHAPEZ=%%(DATA_PATH)s/%%(Z_SAMPLE)sZ_%%(COLOUR)s_galZs.txt\n',
            'DENSZ=%%(DATA_PATH)s/%%(Z_SAMPLE)sZ_galZs.txt\n',
            'SAVE_NZ=%%(DATA_PATH)s/nofz\n',
            '\n',
            '[test]\n',
            'fatal_errors=T\n',
            'save_dir=IAs_%%(TAG)s_%%(Z_SAMPLE)sZ_%%(COLOUR)s\n'
            ]

        os.system('cp %s/wgplus/shell.ini %s/wgplus/%s_like.ini'%(cosmosis_root,cosmosis_root,sample))
        with open("%s/wgplus/%s_like.ini"%(cosmosis_root,sample), "r+") as f:
            old = f.readlines()
	    for line in old:
	    	iniparams.append(line)
            f.seek(0) # rewind
            f.writelines(iniparams)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(fromfile_prefix_chars='@')
    parser.add_argument(
    'datadir',
    help='name of data directory')
    parser.add_argument(
    'tag',
    help='directory specific tag, for saving')
    parser.add_argument(
    'samples',
    nargs='*',
    help="redshift/colour samples for ini-files - 'hzr', 'lzb' etc. or 'all'",
    choices=['all','hzr','hzb','lzr','lzb'])
    parser.add_argument(
    '-param_sampler',
    help="sampler to use - 'test' or 'grid'",
    type=str,
    choices=['test','grid'],
    default='test')
    args = parser.parse_args()

    main(args.datadir,args.tag,args.param_sampler,args.samples)




