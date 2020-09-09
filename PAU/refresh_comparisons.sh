#!/bin/tcsh

foreach colours (LePhare Cigale_2cluster Cigale_3cluster Cigale_3cluster_normi)
	python compare_paus_signals.py ${colours}_signal_comparison \
		OUTPUTS_PAUS_${colours}_uniform_unwindowed \
		OUTPUTS_PAUS_${colours}_uniform_windowed \
		OUTPUTS_PAUS_${colours}_uniform_windowed_qz50 \
		OUTPUTS_PAUS_${colours}_fibonacci_zph-unwindowed \
		OUTPUTS_PAUS_${colours}_fibonacci_zph-windowed
end

foreach pi_rand (uniform_windowed uniform_zph-windowed fibonacci_windowed fibonacci_zph-windowed)
	python compare_paus_signals.py ${pi_rand}_signal_comparison \
		OUTPUTS_PAUS_LePhare_uniform_windowed_qz50 \
		OUTPUTS_PAUS_LePhare_${pi_rand} \
		OUTPUTS_PAUS_Cigale_2cluster_${pi_rand} \
		OUTPUTS_PAUS_Cigale_3cluster_${pi_rand} \
		OUTPUTS_PAUS_Cigale_3cluster_normi_${pi_rand}
end

python compare_paus_signals.py windowed_signal_comparison \
	OUTPUTS_PAUS_LePhare_uniform_windowed_qz50 \
	OUTPUTS_PAUS_LePhare_uniform_windowed OUTPUTS_PAUS_LePhare_fibonacci_zph-windowed \
	OUTPUTS_PAUS_Cigale_2cluster_uniform_windowed OUTPUTS_PAUS_Cigale_2cluster_fibonacci_zph-windowed \
	OUTPUTS_PAUS_Cigale_3cluster_uniform_windowed OUTPUTS_PAUS_Cigale_3cluster_fibonacci_zph-windowed \
	OUTPUTS_PAUS_Cigale_3cluster_normi_uniform_windowed OUTPUTS_PAUS_Cigale_3cluster_normi_fibonacci_zph-windowed

python compare_paus_signals.py unwindowed_signal_comparison \
	OUTPUTS_PAUS_LePhare_uniform_unwindowed_qz50 \
	OUTPUTS_PAUS_LePhare_uniform_unwindowed OUTPUTS_PAUS_LePhare_fibonacci_zph-unwindowed \
	OUTPUTS_PAUS_Cigale_2cluster_uniform_unwindowed OUTPUTS_PAUS_Cigale_2cluster_fibonacci_zph-unwindowed \
	OUTPUTS_PAUS_Cigale_3cluster_uniform_unwindowed OUTPUTS_PAUS_Cigale_3cluster_fibonacci_zph-unwindowed \
	OUTPUTS_PAUS_Cigale_3cluster_normi_uniform_unwindowed OUTPUTS_PAUS_Cigale_3cluster_normi_fibonacci_zph-unwindowed

