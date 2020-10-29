#!/bin/tcsh
#SBATCH -p PREEMPT
#SBATCH --job-name=collect_iterrand
#SBATCH -c 1
#SBATCH --mem=47GB
#SBATCH -t 1-00:00:00

python collect_megacatalogs.py $*
echo ==== OUTFILE: $1
echo ==== FROM:
foreach i (`seq 2 1 $#`)
	echo $argv[$i]
end
echo ====

