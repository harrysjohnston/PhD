#!/bin/tcsh
#SBATCH -p COMPUTE
#SBATCH --job-name=sp_corr
#SBATCH -c 16
#SBATCH --mem=120G
#SBATCH -t 7-00:00:00
source ~/.login
pau

if ($# < 1) then
	echo
	echo 'give job-script directory name'
	echo
	exit 1
endif

set verb = 3
set slop = 0
set base_args = "-verbosity $verb -bin_slop $slop"

set randoms = ('PAUS_KSB_CloneZIDRandoms_6e6_uniqK_Q-0.2.fits'\
			   'PAUS_KSB_CloneZIDRandoms_Unwindowed_uniqK_Q-0.2.fits'\
			   'PAUS_KSB_CloneZIDRandoms_6e6_uniqK_Q-0.2_.pergaldz0.04_zph.fits'\
			   'PAUS_KSB_CloneZIDRandoms_Unwindowed_uniqK_Q-0.2_.pergaldz0.04_zph.fits'\
)
set rand_ids = ('windowed' 'unwindowed' 'zph-windowed' 'zph-unwindowed')
@ Njk = 24
set jobscript_dir = $1

if (-d $jobscript_dir) then
	rm -r $jobscript_dir
endif
mkdir $jobscript_dir
sed "s/JOBSCRIPT_DIR/$jobscript_dir/g" run_jobs.sh > $jobscript_dir/run_jobs.sh

@ x = 0
foreach rand_idx (`seq 1 1 4`)
foreach colors ('LePhare' 'Cigale_2cluster' 'Cigale_3cluster_normi') # 'Cigale_3cluster' 
foreach Pibin ('uniform' 'fibonacci')
#foreach Pibin ('uniform_limrp' 'fibonacci_limrp')
#foreach Pibin ('Pimax140') #('Pimax80' 'Pimax260')# 'finerpbins')
foreach fulldens ('true')# 'false')
foreach qzselect ('false' 'true')
	@ x += 1
	# shorthand colors
	set colors1 = `echo $colors:as/_//:as/igale//:as/luster//:as/hare//:as/e//`
	# only old randoms for uniform-Pi
	#if ($Pibin == 'uniform' && $rand_idx > 2) then
	#	@ x -= 1
	#	continue
	#endif
	# run jackknife 1-by-1 (below)
	set jk_arg = 'jackknife.run=0'
	# set columns
	set r_arg = 'wgplus_config.r_col=comoving_mpc_bcnz_zb'
	# set Pi-binning
	if ($Pibin == 'uniform') then
		set pi_arg = ''
	else if ($Pibin == 'dynamic') then
		set pi_arg = 'wgplus_config.rpar_edges="0.,1.,2.,3.,4.,5.,7.,9.,11.,13.,15.,20.,25.,30.,35.,40.,60.,80.,100.,120.,140.,180.,220."'
	else if ($Pibin == 'fibonacci') then
		set pi_arg = 'wgplus_config.rpar_edges="0.,1.,2.,3.,5.,8.,13.,21.,34.,55.,89.,144.,233."'
	else if ($Pibin == 'fibonacci_finerp') then
		set pi_arg = 'wgplus_config.rpar_edges="0.,1.,2.,3.,5.,8.,13.,21.,34.,55.,89.,144.,233." wgplus_config.max_sep=18. wgplus_config.nbins=7'
	else if ($Pibin == 'uniform_finerp') then
		set pi_arg = 'wgplus_config.max_sep=18. wgplus_config.nbins=7'
	else if ($Pibin == 'fibonacci_coarserp') then
		set pi_arg = 'wgplus_config.rpar_edges="0.,1.,2.,3.,5.,8.,13.,21.,34.,55.,89.,144.,233." wgplus_config.max_sep=18. wgplus_config.nbins=4'
	else if ($Pibin == 'uniform_coarserp') then
		set pi_arg = 'wgplus_config.max_sep=18. wgplus_config.nbins=4'
	else if ($Pibin == 'fibonacci_limrp') then
		set pi_arg = 'wgplus_config.rpar_edges="0.,1.,2.,3.,5.,8.,13.,21.,34.,55.,89.,144.,233." wgplus_config.max_sep=18. wgplus_config.nbins=5'
	else if ($Pibin == 'uniform_limrp') then
		set pi_arg = 'wgplus_config.max_sep=18. wgplus_config.nbins=5'
	else if ($Pibin == 'Pimax80') then
		set pi_arg = 'wgplus_config.min_rpar=-80. wgplus_config.max_rpar=80. wgplus_config.nbin_rpar=40'
	else if ($Pibin == 'Pimax140') then
		set pi_arg = 'wgplus_config.min_rpar=-140. wgplus_config.max_rpar=140. wgplus_config.nbin_rpar=70'
	else if ($Pibin == 'Pimax260') then
		set pi_arg = 'wgplus_config.min_rpar=-260. wgplus_config.max_rpar=260. wgplus_config.nbin_rpar=130'
	endif
	# specify catalogues
	set data_arg = "catalogs.data1=PAUS_KSB.fits"
	set rand_arg = "catalogs.rand1=$randoms[$rand_idx]"
	# select on Qz
	if ($qzselect == 'true') then
		#if ($rand_idx > 2) then
		#	@ x -= 1
		#	continue
		#endif
		set qz1 = '_qz50'
	else if ($qzselect == 'false') then
		set qz1 = ''
	endif
	# choose config file
	set config_r = "PAUS_correlation_template_${colors1}${qz1}_red.ini"
	set config_b = "PAUS_correlation_template_${colors1}${qz1}_blue.ini"
	if ($fulldens == 'true') then
		set config_r = "PAUS_correlation_template_fulldens_${colors1}${qz1}_red.ini"
		set config_b = "PAUS_correlation_template_fulldens_${colors1}${qz1}_blue.ini"
	endif
	# set outdir name
	set rand_id = $rand_ids[$rand_idx]
	set run_id = ${colors}_${Pibin}_${rand_id}${qz1}
	if ($fulldens == 'true') then
		set run_id = ${colors}_fulldens_${Pibin}_${rand_id}${qz1}
	endif
	set outdir = OUTPUTS_PAUS_${run_id}
	# save cats
	if ($x == 1) then
		set save_arg = '-save_cats 1'
	else
		set save_arg = '-save_cats 0'
	endif
	# skip total/wgg correlations
	set idx_args_r = "-index 0 "
	set idx_args_b = "-index 1 "

	set call = "\npython w_pipeline.py $config_r $base_args $idx_args_r $save_arg \
					-p $r_arg $pi_arg $jk_arg $data_arg $rand_arg wgplus_config.save_3d=1 \
					output.savedir=$outdir >& log_PAUS_${run_id}_${x}\n"
	echo $call > $jobscript_dir/array_job_$x.sh
	set call = "\npython w_pipeline.py $config_b $base_args $idx_args_b $save_arg \
					-p $r_arg $pi_arg $jk_arg $data_arg $rand_arg wgplus_config.save_3d=1 \
					output.savedir=$outdir >& log_PAUS_${run_id}_${x}\n"
	echo $call >> $jobscript_dir/array_job_$x.sh
	foreach jkn (`seq 1 1 $Njk`)
		set call = "\npython w_pipeline.py $config_r $base_args $idx_args_r -save_cats 0 \
					-p $r_arg $pi_arg $data_arg $rand_arg jackknife.run=1 jackknife.numbers=$jkn \
					output.savedir=$outdir >& log_PAUS_${run_id}_${x}_jk${jkn}\n"
		echo $call >> $jobscript_dir/array_job_$x.sh
		set call = "\npython w_pipeline.py $config_b $base_args $idx_args_b -save_cats 0 \
					-p $r_arg $pi_arg $data_arg $rand_arg jackknife.run=1 jackknife.numbers=$jkn \
					output.savedir=$outdir >& log_PAUS_${run_id}_${x}_jk${jkn}\n"
		echo $call >> $jobscript_dir/array_job_$x.sh
	end
	chmod +x $jobscript_dir/array_job_$x.sh
end
end
end
end

sbatch -p CORES12 --job-name=PAUSzphcorrs -c 12 --mem=47G -t 7-00:00:00 --array=1-$x%4 $jobscript_dir/run_jobs.sh



