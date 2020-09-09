#!/bin/tcsh

if ( $# < 4 ) then
	echo '\
1: maggies\
2: zmin\
3: zmax\
4: zstep\n'
	exit 1
endif

set maggies = $1
set zmin = $2
set zmax = $3
set zstep = $4
set wd = `pwd`

foreach z ( `seq $zmin $zstep $zmax` )
	echo "#\!/bin/tcsh\
#SBATCH -p COMPUTE\
#SBATCH --job-name=kcorr_z${z}\
#SBATCH -c 1\
#SBATCH --mem=20GB\
#SBATCH -t 0-01:00:00\
source ~/.login\
cd $wd\
\
./shift_maggies.sh $maggies $z\n" > job.sh
	chmod +x job.sh
	ssh splinter-login "cd /share/splinter/hj/PhD/PAU ; sbatch job.sh"
end

