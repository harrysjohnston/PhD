#!/bin/tcsh
#SBATCH -p COMPUTE
#SBATCH --job-name=PAUSrand
#SBATCH -c 16
#SBATCH --mem=120G
#SBATCH --array=1-1%3
#SBATCH -t 3-00:00:00

source ~/.login
pau

set cat = 'PAUS_KSB.fits'
set bl = " "
set zlims = "0.0 1.2"
set dres = 5.
set zph = 1
set zcol = "bcnz_zb"
set diags = 1

set zmax_id = (".pergaldz0.04")

if ($zph == 1) then
	set zmax_id = $zmax_id[$SLURM_ARRAY_TASK_ID]
    set args = "-zph_max $PHD_PATH/PAU/PAUS_KSB${zmax_id}.zmaxtable 320 "
    set Nr = 1
    set suff = '_zph'
else
    set args = " "
    set Nr = 300
    set suff = ''
endif

foreach Q (-0.2)
	foreach kc ('PAUS_uniq.kcorrs')
		set kcid = $kc:as/PAUS_//:as/medians_//:as/.kcorrs//K
		if ($zph == 1) then
			set kcid = "${kcid}_Q${Q}_${zmax_id}${suff}"
		else
			set kcid = "${kcid}_Q${Q}${suff}"
		endif
		set catid = $cat:as/.fits//

		python clone_randoms.py $cat $bl $args -Q $Q -save_diag $diags \
			-magcols mag_i -maglims 22.5 -kcorrs $kc \
			-randoms 1 -idcol numeric_id -zcol $zcol \
			-colnames alpha_j2000 delta_j2000 \
			-id Unwindowed_${kcid} \
			-refresh_zmax 1 -niter 15 -area 19 \
			-zlims $zlims -dres $dres -Nrand $Nr |& grep -v '%' >& log_Unwindowed_${kcid}
		./mask_paus.sh ${catid}_CloneZIDRandoms_Unwindowed_${kcid}.fits
		
		foreach vol (5e6 6e6)
			python clone_randoms.py	$cat $bl $args -Q $Q -save_diag $diags \
				-magcols mag_i -maglims 22.5 -kcorrs $kc \
				-randoms 1 -idcol numeric_id -zcol $zcol \
				-colnames alpha_j2000 delta_j2000 \
				-window $vol -load_windows 0 \
				-id ${vol}_${kcid} \
				-refresh_zmax 0 -niter 15 -area 19 \
				-zlims $zlims -dres $dres -Nrand $Nr |& grep -v '%' >& log_${vol}window_${kcid}
			./mask_paus.sh ${catid}_CloneZIDRandoms_${vol}_${kcid}.fits
		end
	end
end





