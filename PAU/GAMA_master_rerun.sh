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

set savecats = 0
set verb = 3
set slop = 0
set base_args = "-save_cats $savecats -verbosity $verb -bin_slop $slop"

set randoms = ('K1000_SML_CloneZIDRandoms_5e6_gama_McNaughtK_Q-0.7.fits'\
			   'K1000_SML_CloneZIDRandoms_Unwindowed_gama_McNaughtK_Q-0.7.fits'\
			   'K1000_SML_CloneZIDRandoms_5e6_gama_McNaughtK_Q-0.7_pergaldz0.03_zph.fits'\
			   'K1000_SML_CloneZIDRandoms_Unwindowed_gama_McNaughtK_Q-0.7_pergaldz0.03_zph.fits'\
			   '/share/data1/kids/catalogues/randoms/randoms_radec.fits'\
)
set rand_ids = ('windowed' 'unwindowed' 'zph-windowed' 'zph-unwindowed' 'gama-windowed')
set jobscript_dir = $1
@ Njk = 36
set run_jackknife = 3

if (-d $jobscript_dir) then
	rm -r $jobscript_dir
endif
mkdir $jobscript_dir
sed "s/JOBSCRIPT_DIR/$jobscript_dir/g" run_jobs.sh > $jobscript_dir/run_jobs.sh

@ x = 0
foreach rand_idx (`seq 1 1 4`)
foreach zcolname ('zphot_2')# 'Z_TONRY')
foreach Pibin ('uniform' 'fibonacci')# 'uniform_dPi' 'uniform_Pimax' 'uniform_dPi_Pimax')
#foreach Pibin ('uniform_dPi' 'uniform_Pimax' 'uniform_dPi_Pimax')
	@ x += 1
	# short zcolname
	set zcolname1 = `echo $zcolname:as/_// |  tr "[:upper:]" "[:lower:]"`
	# only uniform-Pi for spec-z
	#if ($zcolname == 'Z_TONRY' && $Pibin != 'uniform') then
	#	@ x -= 1
	#	continue
	#endif
	# only GAMA randoms for spec-z & uniform-Pi
	#if ($zcolname == 'Z_TONRY' && $rand_idx != 5) then
	#	@ x -= 1
	#	continue
	#endif
	if (($Pibin == 'uniform_dPi' || $Pibin == 'uniform_Pimax' || $Pibin == 'uniform_dPi_Pimax') && $rand_idx != 5) then
		@ x -= 1
		continue
	endif
	# set columns
	if ($zcolname == 'zphot_2') then
		set r_arg = 'wgplus_config.r_col=comoving_mpc_zphot_2'
	else if ($zcolname == 'Z_TONRY') then
		set r_arg = 'wgplus_config.r_col=comoving_mpc_Z_TONRY'
	endif
	# set Pi-binning
	if ($Pibin == 'uniform') then
		set pi_arg = ''
	else if ($Pibin == 'dynamic') then
		set pi_arg = 'wgplus_config.rpar_edges="0.,1.,2.,3.,4.,5.,7.,9.,11.,13.,15.,20.,25.,30.,35.,40.,60.,80.,100.,120.,140.,180.,220."'
	else if ($Pibin == 'fibonacci') then
		set pi_arg = 'wgplus_config.rpar_edges="0.,1.,2.,3.,5.,8.,13.,21.,34.,55.,89.,144.,233."'
	else if ($Pibin == 'uniform_dPi') then
		set pi_arg = "wgplus_config.nbins_rpar=60"
	else if ($Pibin == 'uniform_Pimax') then
		set pi_arg = "wgplus_config.nbins_rpar=116 wgplus_config.min_rpar=-232. wgplus_config.max_rpar=232."
	else if ($Pibin == 'uniform_dPi_Pimax') then
		set pi_arg = "wgplus_config.nbins_rpar=232 wgplus_config.min_rpar=-232. wgplus_config.max_rpar=232."
	endif
	# specify catalogues
	if ($zcolname == 'Z_TONRY') then
		set data_arg = "catalogs.data1=GAMA_SML_${zcolname1}jk.fits"
	else
		set data_arg = "catalogs.data1=K1000_SML.fits"
	endif
	set rand_arg = "catalogs.rand1=$randoms[$rand_idx]"
	# additional args for GAMA randoms
	if ($rand_idx == 5) then
		set rand_arg = "catalogs.rand1=$randoms[$rand_idx] wgplus_config.rand_r_col=comoving catalogs.rand_cuts1='idmatch(r1, d1, CATAID, CATAID) & (z > 0.02) & (z < 0.5)' catalogs.rand_cuts2='idmatch(r2, d2, CATAID, CATAID) & (z > 0.02) & (z < 0.5)'"
	endif
	# skip total correlations
	set idx_args = "-index 0 1 3 4"
	# make config file
	sed "s/REDSHIFT_COLNAME/$zcolname/g" GAMA_correlation_template.ini > corr_${zcolname}_GAMA_correlations.ini
	# set outdir name
	set rand_id = $rand_ids[$rand_idx]
	set run_id = ${zcolname}_${Pibin}_${rand_id}
	set outdir = OUTPUTS_GAMA_${run_id}

	if ($run_jackknife != 1) then
		# run main correlation
		set call = "\npython w_pipeline.py corr_${zcolname}_GAMA_correlations.ini $base_args $idx_args \
						-p $r_arg $pi_arg $data_arg $rand_arg jackknife.run=0 wgplus_config.save_3d=1 \
						output.savedir=$outdir >& log_GAMA_${run_id}_${x}\n"
		echo $call > $jobscript_dir/array_job_$x.sh
	endif
	# spec-z/other jackknife
	if ($zcolname == 'Z_TONRY' || $run_jackknife != 0) then
		foreach jkn (`seq 1 1 $Njk`)
			set call = "\npython w_pipeline.py corr_${zcolname}_GAMA_correlations.ini $base_args $idx_args \
						-p $r_arg $pi_arg $data_arg $rand_arg jackknife.run=1 jackknife.numbers=$jkn \
						output.savedir=$outdir >& log_GAMA_${run_id}_${x}_jk${jkn}\n"
			echo $call >> $jobscript_dir/array_job_$x.sh
		end
	endif
	chmod +x $jobscript_dir/array_job_$x.sh
end
end
end

sbatch --dependency=afterok:59927 -p CORES16 --job-name=GAMAzphcorrs -c 16 --mem=47G -t 7-00:00:00 --array=1-$x%6 $jobscript_dir/run_jobs.sh



