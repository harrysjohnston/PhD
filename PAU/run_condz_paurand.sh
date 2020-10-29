#!/bin/tcsh
#SBATCH -p COMPUTE
#SBATCH --job-name=CZrand
#SBATCH -c 1
#SBATCH --mem=47G
#SBATCH -t 1-00:00:00

set cat = 'PAUS_KSB_cz.fits'
set Nr = 40
set bl = " "
set zlims = "0.0 1.2"
set args = " "

foreach Q (`seq -1.5 0.1 0.5`)
	foreach kc ('PAUS_uniq.kcorrs')# 'PAUS_medians_5bin.kcorrs' 'PAUS_medians_8bin.kcorrs' 'PAUS_medians_12bin.kcorrs')
		set kcid = $kc:as/PAUS_//:as/medians_//:as/.kcorrs//K
		set kcid = "${kcid}_Q${Q}_cz"

		python clone_randoms.py $cat $bl $args -Q $Q \
			-magcols mag_i -maglims 22.5 -kcorrs $kc \
			-randoms 1 -idcol numeric_id -zcol conditional_z \
			-colnames alpha_j2000 delta_j2000 \
			-id Unwindowed_${kcid} \
			-refresh_zmax 1 -niter 15 -area 19 \
			-zlims $zlims -dres 5. -Nrand $Nr |& grep -v '%' >& log_Unwindowed_${kcid}
		
		foreach vol (4e6)
			python clone_randoms.py	$cat $bl $args -Q $Q \
				-magcols mag_i -maglims 22.5 -kcorrs $kc \
				-randoms 1 -idcol numeric_id -zcol conditional_z \
				-colnames alpha_j2000 delta_j2000 \
				-window $vol -load_windows 0 \
				-id ${vol}_${kcid} \
				-refresh_zmax 0 -niter 15 -area 19 \
				-zlims $zlims -dres 5. -Nrand $Nr |& grep -v '%' >& log_${vol}window_${kcid}
		end
	end
end
