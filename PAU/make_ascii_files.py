# coding: utf-8
from functions import *
if len(sys.argv) < 3:
	print '2x args: flip e1 (1/0) and flip e2 (1/0)'
	sys.exit()

f1 = int(sys.argv[1])
f2 = int(sys.argv[2])
if not (f1 or f2):
	tag = 'f0'
elif f1 and (not f2):
	tag = 'f1'
elif f2 and (not f1):
	tag = 'f2'
elif f1 and f2:
	tag = 'f12'

#dire = 'OUTPUTS_PAUS_KSBtotal_fg1fg2_PW1'
#dire = 'OUTPUTS_PAUS_KSBqzselected_fg1fg2_PW1'
dire = '../OUTPUTS_GAMA_fg1fg2_PW1/'

f1 = [1, -1][int(sys.argv[1])]
f2 = [1, -1][int(sys.argv[2])]

for df in glob('%s/*fits'%dire):
	if 'r2' in df: continue
	print df
	cat = fopen(df)
	if 'PAUS' in dire:
		try:
			ra = np.deg2rad(cat['alpha_j2000'])
			dec = np.deg2rad(cat['delta_j2000'])
			chi = cat['comoving_mpc_bcnz_zb']
			e1 = cat['e1ok'] * f1
			e2 = cat['e2ok'] * f2
			ew = ~np.isnan(e1) & ~np.isnan(e2)
		except:
			ra = np.deg2rad(cat['ra'])
			dec = np.deg2rad(cat['dec'])
			chi = cat['mag_i_cloneComovingDist']
			e1 = np.zeros_like(chi)
			e2 = np.zeros_like(chi)
			ew = np.ones_like(chi)
	elif 'GAMA' in dire:
		try:
			ra = np.deg2rad(cat['RA'])
			dec = np.deg2rad(cat['DEC'])
			chi = cat['comoving_mpc_Z_TONRY']
			e1 = cat['e1_r'] * f1
			e2 = cat['e2_r'] * f2
			ew = ~np.isnan(e1) & ~np.isnan(e2)
		except:
			ra = np.deg2rad(cat['ra'])
			dec = np.deg2rad(cat['dec'])
			chi = cat['PETROMAG_R_cloneComovingDist']
			e1 = np.zeros_like(chi)
			e2 = np.zeros_like(chi)
			ew = np.ones_like(chi)

	tab = np.column_stack((ra, dec, chi, e1, e2, ew))
	if 'qzselect' in dire:
		np.savetxt(df.replace('.fits','_%s.asc'%tag).replace('%s'%dire,'playground/qzselect/'), tab, header='ra\tdec\tchi[mpc/h]\te1\te2\tew')
	elif 'GAMA' in dire:
		np.savetxt(df.replace('.fits','_%s.asc'%tag).replace('%s'%dire,'playground/GAMA/'), tab, header='ra\tdec\tchi[mpc/h]\te1\te2\tew')
	else:
		np.savetxt(df.replace('.fits','_%s.asc'%tag).replace('%s'%dire,'playground'), tab, header='ra\tdec\tchi[mpc/h]\te1\te2\tew')


