#!/bin/tcsh
#SBATCH -p COMPUTE
#SBATCH --job-name=jkclusGAM
#SBATCH --array=0-30%3
#SBATCH -c 16
#SBATCH -t 7-00:00:00
#SBATCH --mem=120GB
source ~/.login
pau

set Cpath = playground/GAMA
set tag = $1
set nbin = $2

python jk_wclust.py \
	$Cpath/gama_red_d1${1}.fits \
	$Cpath/gama_red_d2${1}.fits \
	$Cpath/gama_red_r1${1}.fits \
	$Cpath/gama_red_r2${1}.fits \
	${SLURM_ARRAY_TASK_ID} $1 $2
python jk_wclust.py \
	$Cpath/gama_blue_d1${1}.fits \
	$Cpath/gama_blue_d2${1}.fits \
	$Cpath/gama_blue_r1${1}.fits \
	$Cpath/gama_blue_r2${1}.fits \
	${SLURM_ARRAY_TASK_ID} $1 $2


