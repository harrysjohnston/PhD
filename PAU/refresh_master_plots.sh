#!/bin/tcsh

foreach pausdir (`echo OUTPUTS_PAUS_*uniform* OUTPUTS_PAUS_*fibonacci*`)
foreach gamadir (`echo OUTPUTS_GAMA_*uniform* OUTPUTS_GAMA_*fibonacci*`)
	#if ($pausdir =~ *qz50*) then
		#continue
	python plot_master_signal_compar.py -gama OUTPUTS_GAMA_ztonry_old $gamadir -paus $pausdir
		#python plot_master_signal_compar.py -gama OUTPUTS_GAMA_ztonry_old OUTPUTS_GAMA_zphot_2_fibonacci_zph-windowed -paus $pausdir
	#endif
end
end

