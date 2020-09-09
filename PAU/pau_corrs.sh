#!/bin/tcsh
#SBATCH -p SMP
#SBATCH --job-name=paucorr
#SBATCH -c 40
#SBATCH --mem=120G
#SBATCH -t 7-00:00:00
source ~/.login
pau

set config = pau_total_config.ini
set confi2 = pau_qzselected_config.ini
set confi3 = pau_vmaxweighted_config.ini
set confi4 = pau_vmax_qz_config.ini
set confi5 = pau_lowMr_config.ini
set confi6 = pau_midMr_config.ini
set confi7 = pau_highMr_config.ini
set confi8 = pau_area_config.ini
set confi9 = pau_lowz_config.ini
set confi10 = pau_qzweighted_config.ini
set bs = 0
set vs = 2
set sc = 0
set save_3d = 0
# 1=jk-only, 2=jk-first, 3=jk-last, 4=jk-collect
set jk = 1
set estim = PW1
set s1 = "jackknife.run=$jk\
wgplus_config.save_3d=$save_3d\
wgplus_config.min_sep=0.1\
wgplus_config.max_sep=50.\
wgplus_config.nbins=6\
wgplus_config.min_rpar=-60.\
wgplus_config.max_rpar=60.\
wgplus_config.nbins_rpar=30\
wgplus_config.estimator=$estim\
wgplus_config.compensated=1"
set fg1 = wgplus_config.flip_g1
set fg2 = wgplus_config.flip_g2
set args = "-save_cats $sc -verbosity $vs -bin_slop $bs"

# skip clustering for e1/2 flips
foreach conf ($config)
	set id1 = $conf:as/pau_//:as/randomstest//:as/_config.ini//
	set id = "_fg1fg2_${estim}"
	python w_pipeline.py $conf -num_threads 40 -index 2 4 $args -p `echo $s1` "$fg1"=1 "$fg2"=1 output.savedir="OUTPUTS_PAUS_KSB${id1}${id}" >& logwpipe_`ts`_PAUS_KSB${id1}${id}
	#set id = "_f0_${estim}"
	#python w_pipeline.py $conf -index 0 2 4 $args -p `echo $s1` "$fg1"=0 "$fg2"=0 output.savedir="OUTPUTS_PAUS_KSB${id1}${id}" >& logwpipe_`ts`_PAUS_KSB${id1}${id}
	#set id = "_fg2_${estim}"
	#python w_pipeline.py $conf -index 0 2 4 $args -p `echo $s1` "$fg1"=1 "$fg2"=0 output.savedir="OUTPUTS_PAUS_KSB${id1}${id}" >& logwpipe_`ts`_PAUS_KSB${id1}${id}
	#set id = "_fg1_${estim}"
	#python w_pipeline.py $conf -index 0 2 4 $args -p `echo $s1` "$fg1"=0 "$fg2"=1 output.savedir="OUTPUTS_PAUS_KSB${id1}${id}" >& logwpipe_`ts`_PAUS_KSB${id1}${id}
end


