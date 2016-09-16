from __future__ import print_function, division
from os.path import join, normpath, basename, isdir
import os
import argparse

def main(rand_cat, qsub_dir):
	
	if not isdir(qsub_dir):
		os.mkdir(qsub_dir)

	for i in range(1,11):
		shell_script = [
		'#!/bin/tcsh',
		'#PBS -q compute',
		'#PBS -N fits2ascii',
		'#PBS -l nodes=1:ppn=1',
		'#PBS -l mem=30gb',
		'#PBS -l walltime=120:00:00',
		'#PBS -M zcaphjo@ucl.ac.uk',
		'#PBS -m a',
		'',
		'module load dev_tools/nov2014/python-anaconda',
		'',
		'cd $PBS_O_WORKDIR',
		'',
		'python /share/splinter/ug_hj/PhD/fits2ascii.py ' + str(rand_cat) + ' @/share/splinter/ug_hj/PhD/random_cols /share/splinter/ug_hj/PhD/steps_randcats --step ' + str(i)
		]

		shell_script.append("")

		F = join(qsub_dir, 'step' + str(i).zfill(2) + "_qsub.sh")
		A = open(F, "w")
		T = "\n".join(shell_script)
		A.write(str(T))
		A.close()

			# os.system("qsub " + F)
			# os.system("sleep 1")

if __name__ == "__main__":
	parser = argparse.ArgumentParser()
	parser.add_argument(
		'rand_cat',
		help='random fits catalog to be converted to ascii in parallel steps')
	parser.add_argument(
		'qsub_dir',
		help='destination directory for generated shell scripts')
	args = parser.parse_args()
	main(args.rand_cat, args.qsub_dir)
