#############################################################
### RUN FROM A DIRECTORY WHERE YOU WANT THE RESULTS FILES ###
#############################################################
import os
from os.path import join

keys = ['z2_b','z2_r','z1_b','z1_r','sdss_b','sdss_r']
file_names = dict(zip(keys, ['highZ_Blue', 'highZ_Red', 'lowZ_Blue', 'lowZ_Red', 'highZ_Blue', 'highZ_Red']))

# SPECIFY PATHS TO POINT TO CATALOG_SAMPLER OUTPUT DIRECTORY
PHD = "/share/splinter/hj/PhD/"
# Standard
#gama = PHD + "Wcorr_jk3x5_24_02_18_z0.22_c0.66/"
gama = PHD + "Wcorr_jk3x5_24_02_18_z0.26_c0.66/"
sdss = PHD + "Wcorr_SDSS_jk10x10_25_02_18_c0.66/"
# HH
HH = 0
if HH:
	gama = PHD + "Wcorr_HH_jk4x5_25_02_18_allz/"
	sdss = PHD + "Wcorr_HH_SDSS_jk10x10_25_02_18_allz/"
# largePi
largePi = 1
if largePi:
	gama = PHD + "Wcorr_GAMA_largePi_jk3x5_6rpbins_z0.26_c0.66/"
	sdss = PHD + "Wcorr_SDSS_largePi_jk10x10_6rpbins_c0.66/"

# COPY WG+/x FILES FROM GAMA AND/OR SDSS OUTPUT DIRS
do_gama = 1
do_sdss = 1

IA = 1
clus = 0

wgx = 0
covar_pref = ['JKcovarP_', 'JKcovarX_'][wgx]
if wgx:
	IA = 1
	largePi = 0
	clus = 0

dirs = dict(zip(keys, [gama]*4 + [sdss]*2))
clus_prefs = dict(zip(keys, ['MICE']*4 + ['swot']*2))
files = {}
for k in keys:
	files['wgp_'+k] = file_names[k].split('_')[0] + '_vs_' + file_names[k]
	files['covar_wgp_'+k] = covar_pref + file_names[k]
	if largePi:
		files['wgp_'+k] += '_largePi'
		files['covar_wgp_'+k] += '_largePi'

	files['wgg_'+k] = 'wgg_%s_'%clus_prefs[k] + file_names[k] + '_UnMasked'
	files['covar_wgg_'+k] = 'wggJK_covar_%s_'%clus_prefs[k] + file_names[k] + '_UnMasked'

if do_gama: print(gama)
if do_sdss: print(sdss)
if do_gama & (not do_sdss): keys = keys[:4]
if do_sdss & (not do_gama): keys = keys[-2:]

for k in keys:
	if HH & (k not in ['z2_r', 'sdss_r']):
		continue
	if IA:
		#print("%s -> wgp_%s" % (files['wgp_'+k], k))
		os.system("cp %s ./wgp_%s" % (join(dirs[k], 'to_plot', files['wgp_'+k]), k))

		#print("%s -> covar_wgp_%s" % (files['covar_wgp_'+k], k))
		os.system("cp %s ./covar_wgp_%s" % (join(dirs[k], 'to_plot', files['covar_wgp_'+k]), k))
	if clus:
		#print("%s -> wgg_%s" % (files['wgg_'+k], k))
		os.system("cp %s ./wgg_%s" % (join(dirs[k], 'to_plot', files['wgg_'+k]), k))

		#print("%s -> covar_wgg_%s" % (files['covar_wgg_'+k], k))
		os.system("cp %s ./covar_wgg_%s" % (join(dirs[k], 'to_plot', files['covar_wgg_'+k]), k))
	


