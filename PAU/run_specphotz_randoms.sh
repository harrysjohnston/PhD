#!/bin/tcsh
#SBATCH -p PREEMPT
#SBATCH --job-name=sp_rand
#SBATCH -c 1
#SBATCH --mem=47G
#SBATCH -t 3-00:00:00

source ~/.login
pau

set cat = "$PHD_PATH/PAU/K1000_SML.fits"
set catid = $cat:as/.fits//
set args = " "
set zlims = "0.0 0.6"
set Nr = 50
set dres = 5.
set Q = '-0.7'
#set id_ph = "phot_rand"

set zph = 1

if ($zph == 1) then
	set args = "-zph_max $PHD_PATH/PAU/K1000_SML.zmaxtable 100 "
	set Nr = 4
	set suff = '_zph'
else
	set args = " "
	set Nr = 40
	set suff = ''
endif

# CHECK KCORRS!
set kc = "$PHD_PATH/SMLambdarApMatchedPhotom.zphot2.kcorrs_h5"
set id_sp = "spec_rand"
set id_ph1 = "zphot_1"
set id_ph2 = "zphot_2"

foreach id ($id_ph2)# $id_ph2)
	set zid = $id
	set id = ${id}_uniqk
	python clone_randoms.py $cat $args -Q $Q \
		-magcols PETROMAG_R -maglims 19.8 -kcorrs $kc \
		-randoms 1 -idcol CATAID -zcol $zid \
		-colnames RA DEC \
		-id Unwindowed_${id} \
		-refresh_zmax 1 -niter 15 -area 180. \
		-zlims $zlims -dres $dres -Nrand $Nr |& grep -v '%' >& log_Unwindowed_${id}
	python $PHD_PATH/get_gama_rand_radec.py ${catid}_CloneZIDRandoms_Unwindowed_${id}.fits

	foreach vol (4e6)
		python clone_randoms.py $cat $args -Q $Q \
			-magcols PETROMAG_R -maglims 19.8 -kcorrs $kc \
			-randoms 1 -idcol CATAID -zcol $zid \
			-colnames RA DEC \
			-window ${vol} -load_windows 0 \
			-id ${vol}_${id} \
			-refresh_zmax 1 -niter 15 -area 180. \
			-zlims $zlims -dres $dres -Nrand $Nr |& grep -v '%' >& log_${vol}_${id}
		python $PHD_PATH/get_gama_rand_radec.py ${catid}_CloneZIDRandoms_${vol}_${id}.fits
	end
end

#python clone_randoms.py $cat $args -Q $Q \
#	-magcols PETROMAG_R -maglims 19.8 -kcorrs $kc \
#	-randoms 1 -idcol CATAID -zcol Z_TONRY \
#	-colnames RA DEC \
#	-id Unwindowed_${id_sp} \
#	-refresh_zmax 1 -niter 15 -area 180. \
#	-zlims $zlims -dres $dres -Nrand $Nr |& grep -v '%' >& log_Unwindowed_${id_sp}
#python $PHD_PATH/get_gama_rand_radec.py ${catid}_CloneZIDRandoms_Unwindowed_${id_sp}.fits
#
#set vol = 4e6
#python clone_randoms.py $cat $args -Q $Q \
#	-magcols PETROMAG_R -maglims 19.8 -kcorrs $kc \
#	-randoms 1 -idcol CATAID -zcol Z_TONRY \
#	-colnames RA DEC \
#	-id ${vol}_${id_sp} \
#	-window $vol -load_windows 0 \
#	-refresh_zmax 1 -niter 15 -area 180. \
#	-zlims $zlims -dres $dres -Nrand $Nr |& grep -v '%' >& log_${vol}_${id_sp}
#python $PHD_PATH/get_gama_rand_radec.py ${catid}_CloneZIDRandoms_${vol}_${id_sp}.fits
#
#python clone_randoms.py $cat $args -Q $Q \
#	-magcols PETROMAG_R -maglims 19.8 -kcorrs $kc \
#	-randoms 1 -idcol CATAID -zcol z_ANNZ_KV \
#	-colnames RA DEC \
#	-id Unwindowed_${id_ph} \
#	-refresh_zmax 1 -niter 15 -area 180. \
#	-zlims $zlims -dres $dres -Nrand $Nr |& grep -v '%' >& log_Unwindowed_${id_ph}
#python $PHD_PATH/get_gama_rand_radec.py ${catid}_CloneZIDRandoms_Unwindowed_${id_ph}.fits

# no k-corrections to make this work.....
#set cat = "$PHD_PATH/PAU/PAUS_KSB.fits"
#set catid = $cat:as/.fits//
##set args = " "
#set zlims = "0.0 1.2"
##set Nr = 50
##set dres = 5.
##set zph = 1
#set Q = '-0.2'
#set kc = "$PHD_PATH/PAU/PAUS_uniq.kcorrs"
#
#foreach id ('pz_draw' 'sz_draw')
#	python clone_randoms.py $cat $args -Q $Q \
#		-magcols mag_i -maglims 22.5 -kcorrs $kc \
#		-randoms 1 -idcol numeric_id -zcol $id \
#		-colnames alpha_j2000 delta_j2000 \
#		-id Unwindowed_${id} \
#		-refresh_zmax 1 -niter 15 -area 19. \
#		-zlims $zlims -dres $dres -Nrand $Nr |& grep -v '%' >& log_Unwindowed_${id}
#	./mask_paus.sh PAUS_KSB_CloneZIDRandoms_Unwindowed_${id}.fits
#end


