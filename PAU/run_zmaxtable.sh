#!/bin/tcsh
#SBATCH -p COMPUTE
#SBATCH --job-name=zmax
#SBATCH --reservation=hj_3
#SBATCH -c 16
#SBATCH --array=2-2%2
#SBATCH --mem=120G
#SBATCH -t 7-00:00:00

set Ndraw = 640
set gdz = 0.03
set pdz = 0.04
set mode = 'pergal'

if ($SLURM_ARRAY_TASK_ID == 1) then
	python make_zmax_table.py K1000_SML.fits zphot_2 Z_TONRY CATAID $Ndraw ../SMLambdarApMatchedPhotom_McNaught.kcorrs 19.8 PETROMAG_R '-0.7' 0.6 $gdz $mode >& log_zmaxtable_gama_dz${gdz}
else if ($SLURM_ARRAY_TASK_ID == 2) then
	python make_zmax_table.py PAUS_KSB.fits bcnz_zb ZBEST numeric_id $Ndraw PAUS_uniq.kcorrs 22.5 mag_i '-0.2' 1.2 $pdz $mode >& log_zmaxtable_paus_dz${pdz}
endif



