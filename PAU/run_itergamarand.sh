#!/bin/tcsh
#SBATCH -p PREEMPT
#SBATCH --job-name=GAMArand
#SBATCH -c 1
#SBATCH --array=1-320%16
#SBATCH --mem=20GB
#SBATCH -t 3-00:00:00

source ~/.login
pau

set mcatid = gamacat
python make_megacatalogs.py K1000_SML.pergaldz0.03.zmaxtable K1000_SML.fits zphot_2 PETROMAG_R CATAID RA DEC $SLURM_ARRAY_TASK_ID $mcatid

set ind = `printf "%03d\n" ${SLURM_ARRAY_TASK_ID}`
set cat = ${mcatid}${ind}.fits
set bl = " "
set zlims = "0.0 0.6"
set dres = 5.
set zph = 0
set zcol = "zspec"
set zmaxcol = "zmax"

set args = "-zmax_col zmax "
set Nr = 1
set suff = ''
#endif

foreach Q (-0.7)
	foreach kc ("$PHD_PATH/SMLambdarApMatchedPhotom_McNaught.kcorrs")
		set kcid = 'gama_McNaughtK'
		if ($zph == 1) then
			set kcid = "${kcid}_Q${Q}_${nbin_zph}${suff}"
		else
			set kcid = "${kcid}_Q${Q}${suff}"
		endif
		set catid = $cat:as/.fits//

		python clone_randoms.py $cat $bl $args -Q $Q \
			-magcols PETROMAG_R -maglims 19.8 -kcorrs $kc \
			-randoms 1 -idcol CATAID -zcol $zcol \
			-colnames RA DEC \
			-id Unwindowed_${kcid}_${ind} \
			-refresh_zmax 1 -niter 15 -area 180. \
			-zlims $zlims -dres $dres -Nrand $Nr |& grep -v '%' >& log_Unwindowed_${kcid}_${ind}
		python $PHD_PATH/get_gama_rand_radec.py ${catid}_CloneZIDRandoms_Unwindowed_${kcid}.fits

		foreach vol (4e6)
			python clone_randoms.py	$cat $bl $args -Q $Q \
				-magcols PETROMAG_R -maglims 19.8 -kcorrs $kc \
				-randoms 1 -idcol CATAID -zcol $zcol \
				-colnames RA DEC \
				-window $vol -load_windows 0 \
				-id ${vol}_${kcid}_${ind} \
				-refresh_zmax 0 -niter 15 -area 180. \
				-zlims $zlims -dres $dres -Nrand $Nr |& grep -v '%' >& log_${vol}window_${kcid}_${ind}
			python $PHD_PATH/get_gama_rand_radec.py ${catid}_CloneZIDRandoms_${vol}_${kcid}.fits
		end
	end
end
rm $cat



