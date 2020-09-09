#!/bin/tcsh
#SBATCH -p COMPUTE
#SBATCH --job-name=kcorr_z0.650
#SBATCH -c 1
#SBATCH --mem=20GB
#SBATCH -t 0-01:00:00
source ~/.login
cd /share/splinter/hj/PhD/PAU

./shift_maggies.sh /share/splinter/hj/PhD/SMLambdarApMatchedPhotom.maggies.dat 0.650

