#!/bin/tcsh
#SBATCH -p PREEMPT
#SBATCH --job-name=PAUSrand
#SBATCH -c 1
#SBATCH --array=1-320%16
#SBATCH --mem=20GB
#SBATCH -t 3-00:00:00

source ~/.login
pau

set mcatid = paucat
python make_megacatalogs.py PAUS_KSB.pergaldz0.04.zmaxtable PAUS_KSB.fits bcnz_zb mag_i numeric_id alpha_j2000 delta_j2000 $SLURM_ARRAY_TASK_ID $mcatid

set ind = `printf "%03d\n" ${SLURM_ARRAY_TASK_ID}`
set cat = ${mcatid}${ind}.fits
set bl = " "
set zlims = "0.0 1.2"
set dres = 5.
set zph = 0
set zcol = "zspec"
set zmaxcol = "zmax"

set args = "-zmax_col zmax "
set Nr = 1
set suff = ''
#endif

foreach Q (-0.2)
	foreach kc ('PAUS_uniq.kcorrs')
		set kcid = $kc:as/PAUS_//:as/medians_//:as/.kcorrs//K
		if ($zph == 1) then
			set kcid = "${kcid}_Q${Q}_${nbin_zph}${suff}"
		else
			set kcid = "${kcid}_Q${Q}${suff}"
		endif
		set catid = $cat:as/.fits//

		python clone_randoms.py $cat $bl $args -Q $Q \
			-magcols mag_i -maglims 22.5 -kcorrs $kc \
			-randoms 1 -idcol numeric_id -zcol $zcol \
			-colnames alpha_j2000 delta_j2000 \
			-id Unwindowed_${kcid}_${ind} \
			-refresh_zmax 1 -niter 15 -area 19 \
			-zlims $zlims -dres $dres -Nrand $Nr |& grep -v '%' >& log_Unwindowed_${kcid}_${ind}
		./mask_paus.sh ${catid}_CloneZIDRandoms_Unwindowed_${kcid}.fits

		foreach vol (5e6)
			python clone_randoms.py	$cat $bl $args -Q $Q \
				-magcols mag_i -maglims 22.5 -kcorrs $kc \
				-randoms 1 -idcol numeric_id -zcol $zcol \
				-colnames alpha_j2000 delta_j2000 \
				-window $vol -load_windows 0 \
				-id ${vol}_${kcid}_${ind} \
				-refresh_zmax 0 -niter 15 -area 19 \
				-zlims $zlims -dres $dres -Nrand $Nr |& grep -v '%' >& log_${vol}window_${kcid}_${ind}
			./mask_paus.sh ${catid}_CloneZIDRandoms_${vol}_${kcid}.fits
		end
	end
end
rm $cat



