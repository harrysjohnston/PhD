from __future__ import print_function,division
import numpy as np
from os.path import join
import os

def main(datadir,tag,param_sampler,*samples):
    cosmosis_root = '/share/splinter/hj/PhD/cosmosis'
    if samples=='all':
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
            '[runtime]',
            'sampler = grid',
            '',
            '[DEFAULT]',
            'N_BIN=50',
            'N_Z=50',
            'MY_ROOT=/share/splinter/hj/PhD',
            'DATA_PATH=%%(MY_ROOT)s/%s'%datadir,
            'TAG=%s'%tag,
            'Z_SAMPLE=%s'z,
            'COLOUR=%s'col,
            'SHAPEZ=%%(DATA_PATH)s/%%(Z_SAMPLE)sZ_%%(COLOUR)s_galZs.txt',
            'DENSZ=%%(DATA_PATH)s/%%(Z_SAMPLE)sZ_galZs.txt',
            'SAVE_NZ=%%(DATA_PATH)s/nofz',
            '',
            '[grid]',
            'fatal_errors=T',
            'nsample_dimension=10',

            '[output]',
            'format=text',
            'filename=IAs_%%(TAG)s_likelihoodtest'
            ]
        elif param_sampler=='test':
            iniparams = [
            '[runtime]',
            'sampler = test',
            '',
            '[DEFAULT]',
            'N_BIN=50',
            'N_Z=50',
            'MY_ROOT=/share/splinter/hj/PhD',
            'DATA_PATH=%%(MY_ROOT)s/%s'%datadir,
            'TAG=%s'%tag,
            'Z_SAMPLE=%s'z,
            'COLOUR=%s'col,
            'SHAPEZ=%%(DATA_PATH)s/%%(Z_SAMPLE)sZ_%%(COLOUR)s_galZs.txt',
            'DENSZ=%%(DATA_PATH)s/%%(Z_SAMPLE)sZ_galZs.txt',
            'SAVE_NZ=%%(DATA_PATH)s/nofz',
            '',
            '[test]',
            'fatal_errors=T',
            'save_dir=IAs_%%(TAG)s'
            ]

        os.system('cp %s/wgplus/shell.ini %s/wgplus/%s_like.ini'%(cosmosis_root,cosmosis_root,sample))
        with open("%s/wgplus/%s_like.ini"%(cosmosis_root,sample), "r+") as f:
            old = f.read() # read everything in the file
            f.seek(0) # rewind
            f.write(iniparams + old)

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




