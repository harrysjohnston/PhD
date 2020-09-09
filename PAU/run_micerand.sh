#!/bin/tcsh
#SBATCH -p COMPUTE
#SBATCH --job-name=MICErand
#SBATCH -c 1
#SBATCH --mem=13G
#SBATCH -t 1-00:00:00

set cat = 'MICE_forPAUSrandomstesting_dec20.fits'
set Nr = 40
set bl = " "
set args = " "

foreach Q (`seq -1.5 0.1 1.5`)
	foreach kc ('MICE_dec20_4bin.kcorrs')
		set kcid = $kc:as/MICE_//:as/dec20_//:as/.kcorrs//K
		set kcid = "${kcid}_Q${Q}"

		python clone_randoms.py $cat $bl $args -Q $Q \
			-magcols sdss_i_true -maglims 22.5 -kcorrs $kc -randoms 1 \
			-colnames ra_gal dec_gal -area 180 -zcol z_cgal -idcol unique_gal_id \
			-dres 5 -zlims 0.0 1.2 -niter 15 -Nrand $Nr \
			-id Unwindowed_${kcid} -refresh_zmax 1 |& grep -v '%' >& log_miceUnwindowed_${kcid}

		foreach vol (4e6 6e6)
			python clone_randoms.py $cat $bl $args -Q $Q \
				-magcols sdss_i_true -maglims 22.5 -kcorrs $kc -randoms 1 \
				-colnames ra_gal dec_gal -area 180 -zcol z_cgal -idcol unique_gal_id \
				-dres 5 -zlims 0.0 1.2 -niter 15 -Nrand $Nr \
				-window ${vol} -id mice${vol}_${kcid} -load_windows 0 -refresh_zmax 0 |& grep -v '%' >& log_mice${vol}_${kcid}
		end
	end
end
