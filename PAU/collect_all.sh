#!/bin/tcsh

foreach x ('Unwindowed' '4e6')

	set outcat = gamacat_CloneZIDRandoms_${x}_gama_McNaughtK_Q-0.7.fits
	set incats = `echo "gamacat*_CloneZIDRandoms_${x}_gama_McNaughtK_Q-0.7_*.fits"`

	sbatch collect_iterrand.sh $outcat "$incats"
	#echo collect_iterrand.sh $outcat "$incats"
end

foreach x ('Unwindowed' '5e6')

	set outcat = paucat_CloneZIDRandoms_${x}_uniqK_Q-0.2.fits
	set incats = `echo "paucat*_CloneZIDRandoms_${x}_uniqK_Q-0.2_*.fits"`

	sbatch collect_iterrand.sh $outcat "$incats"
	#echo collect_iterrand.sh $outcat "$incats"
end

