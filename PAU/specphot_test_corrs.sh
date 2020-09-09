#!/bin/tcsh
#SBATCH -p COMPUTE
#SBATCH --job-name=sp_corr
#SBATCH -c 16
#SBATCH --array=1-1%1
#SBATCH --mem=120G
#SBATCH -t 7-00:00:00
source ~/.login
pau

set savecats = 0
set verb = 1
set slop = 0.01
set args = "-save_cats $savecats -verbosity $verb -bin_slop $slop"
# 1=jk-only, 2=jk-first, 3=jk-last, 4=jk-collect
set s1 = "jackknife.run=3\
wgplus_config.save_3d=1\
wgplus_config.min_sep=0.1\
wgplus_config.max_sep=60.\
wgplus_config.nbins=11\
wgplus_config.min_rpar=-60.\
wgplus_config.max_rpar=60.\
wgplus_config.nbins_rpar=30\
wgplus_config.compensated=1"

if ($# < 6) then
	echo "\ngive:\
1 = mode: 1 for clustering, 2 for IA \
2 = path to outdir \
3 = windowed (give vol string e.g. 4e6) or unwindowed (0) \
4 = uniform Pi-binning (0) or dynamic (1)\
5 = McNaught k-corrs (0) or unique (1)\
6 = use standard method (0) or drawn z-spec method (1)"
	echo
	exit
endif

set mode = $1

if ($6 == 0) then
	set confs = (DspecRspec_config.ini Dzphot_2Rspec_config.ini Dzphot_2Rzphot_2_config.ini Dzphot_1Rspec_config.ini)
else if ($6 == 1) then
	set confs = (Dzphot_2Rzphot_2_zph_config.ini)
	if ($SLURM_ARRAY_TASK_ID > 1) then
		echo "not re-running Rspec correlations for drawn z-spec methods"
		exit
	endif
endif

set conf = $confs[${SLURM_ARRAY_TASK_ID}]
#set conf = $confs[1]
##python skyknife.py $conf >& logskyknife_`ts`_${id}
set id = $conf:as/_config.ini//

if ($3 == 0) then
	set c1 = corr_$conf
	cp $conf $c1
else
	set c1 = corr_win_$conf
	sed "s/Unwindowed/$3/g" $conf > $c1
endif

if ($5 == 1) then
	sed -i 's/spec_rand.fits/spec_rand_uniqk.fits/g' $c1
	sed -i 's/zphot_2.fits/zphot_2_uniqk.fits/g' $c1
	sed -i 's/zphot_1.fits/zphot_1_uniqk.fits/g' $c1
endif

if ($4 == 0) then
	set arg = ''
else if ($4 == 1) then
	set arg = 'wgplus_config.rpar_edges="0.,1.,2.,3.,4.,5.,7.,9.,11.,13.,15.,20.,25.,30.,35.,40.,60.,80.,100.,120.,140.,180.,220."'
	#set arg = 'wgplus_config.rpar_edges="0.,1.,2.,3.,5.,8.,13.,21.,34.,55.,89.,144.,233."' # fibonacci
endif

if ($mode == 1) then
	python w_pipeline.py $c1 $args \
		-p `echo $s1` output.savedir=$2 $arg \
		-index 0 1 2 >& logwpipe_`ts`_GAMA_${id}
endif

if ($mode == 2) then
	python w_pipeline.py $c1 $args \
		-p `echo $s1` output.savedir=$2 $arg \
		-index 3 4 5 >& logwpipe_`ts`_GAMA_${id}
	python w_pipeline.py $c1 -p jackknife.run=4
	#cp OUTPUTS_specphot_test/wgp_${id}*dat OUTPUTS_specphot_test/f0/
endif

	#python w_pipeline.py $c1 $args -p `echo $s1` wgplus_config.flip_g1=1 wgplus_config.flip_g2=1 -index 0 >& logwpipe_`ts`_GAMA_${id}
	##python w_pipeline.py $c1 -p jackknife.run=4
	#cp OUTPUTS_specphot_test/wgp_${id}*dat OUTPUTS_specphot_test/f12/

	#python w_pipeline.py $c1 $args -p `echo $s1` wgplus_config.flip_g1=1 wgplus_config.flip_g2=0 -index 0 >& logwpipe_`ts`_GAMA_${id}
	##python w_pipeline.py $c1 -p jackknife.run=4
	#cp OUTPUTS_specphot_test/wgp_${id}*dat OUTPUTS_specphot_test/f1/

	#python w_pipeline.py $c1 $args -p `echo $s1` wgplus_config.flip_g1=0 wgplus_config.flip_g2=1 -index 0 >& logwpipe_`ts`_GAMA_${id}
	##python w_pipeline.py $c1 -p jackknife.run=4
	#cp OUTPUTS_specphot_test/wgp_${id}*dat OUTPUTS_specphot_test/f2/
