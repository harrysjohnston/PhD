#!/bin/tcsh
#SBATCH -p COMPUTE
#SBATCH --job-name=GAMArand
#SBATCH -c 16
#SBATCH --mem=120G
#SBATCH --array=1-1%1
#SBATCH -t 3-00:00:00

source ~/.login
pau

set cat = "$PHD_PATH/PAU/K1000_SML.fits"
set bl = " "
set zlims = "0.0 0.6"
set dres = 5.
set zph = 1
set zcol = "zphot_2"

set zmax_id = ("pergaldz0.03")# "pergaldz0.005")

if ($zph == 1) then
	set zmax_id = $zmax_id[$SLURM_ARRAY_TASK_ID]
	set args = "-zph_max $PHD_PATH/PAU/K1000_SML.${zmax_id}.zmaxtable 320 "
	set Nr = 1
	set suff = '_zph'
else
	set args = " "
	set Nr = 300
	set suff = ''
endif

foreach Q (-0.7)
	foreach kc ("$PHD_PATH/SMLambdarApMatchedPhotom_McNaught.kcorrs")
		set kcid = 'gama_McNaughtK'
		if ($zph == 1) then
			set kcid = "${kcid}_Q${Q}_${zmax_id}${suff}"
		else
			set kcid = "${kcid}_Q${Q}${suff}"
		endif
		set catid = $cat:as/.fits//

		python clone_randoms.py $cat $bl $args -Q $Q \
			-magcols PETROMAG_R -maglims 19.8 -kcorrs $kc \
			-randoms 1 -idcol CATAID -zcol $zcol \
			-colnames RA DEC \
			-id Unwindowed_${kcid} \
			-refresh_zmax 1 -niter 15 -area 180. \
			-zlims $zlims -dres $dres -Nrand $Nr |& grep -v '%' >& log_Unwindowed_${kcid}
		python $PHD_PATH/get_gama_rand_radec.py ${catid}_CloneZIDRandoms_Unwindowed_${kcid}.fits

		foreach vol (5e6)
			python clone_randoms.py	$cat $bl $args -Q $Q \
				-magcols PETROMAG_R -maglims 19.8 -kcorrs $kc \
				-randoms 1 -idcol CATAID -zcol $zcol \
				-colnames RA DEC \
				-window $vol -load_windows 1 \
				-id ${vol}_${kcid} \
				-refresh_zmax 0 -niter 15 -area 180. \
				-zlims $zlims -dres $dres -Nrand $Nr |& grep -v '%' >& log_${vol}window_${kcid}
			python $PHD_PATH/get_gama_rand_radec.py ${catid}_CloneZIDRandoms_${vol}_${kcid}.fits
		end
	end
end

