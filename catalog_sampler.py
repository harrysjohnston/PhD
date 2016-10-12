from __future__ import print_function, division
from astropy.io import fits
import numpy as np
from astropy.io import ascii
from os.path import join, isdir, basename, normpath
from os import listdir, mkdir
import os
import argparse
import csv
from astropy import cosmology
from astropy.cosmology import Planck13
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import sys
import scipy.integrate as scint
import scipy.stats as stat
from scipy.stats import chi2

class RealCatalogue:

	def __init__(self, path):
		"""""
		read-in catalogue

		"""""
		self.path = path
		hdulist = fits.open(path)
		self.data = hdulist[1].data
		self.columns = hdulist[1].columns.names
		# ascii (.asc) file IDs
		self.labels = ['highZ_Red', 'highZ_Blue', 'lowZ_Red', 'lowZ_Blue', 'highZ', 'lowZ']
		# bodies of wcorr output file IDs
		self.wcorrLabels = ['hiZ_vs_hiZ_R', 'hiZ_vs_hiZ_B', 'loZ_vs_loZ_R', 'loZ_vs_loZ_B']

	def cut_data(self, pgm_, z_, colour_, *bitmask_): 
		"""""
		cut catalogue according to bitmasks, 
		PGM, & into subsamples
		
		"""""
		assert 'RA_1' in self.columns, "'RA_1' not in columns, see column headers: "+ str(self.columns)
		assert 'DEC_1' in self.columns, "'DEC_1' not in columns, see column headers: "+ str(self.columns)
		assert 'pgm' in self.columns, "'pgm' not in columns, see column headers: "+ str(self.columns)
		assert 'Z_1_1' in self.columns, "'Z_1_1' not in columns, see column headers: "+ str(self.columns)
		assert 'absmag_g_1' in self.columns, "'absmag_g_1' not in columns, see column headers: "+ str(self.columns)
		assert 'absmag_i_1' in self.columns, "'absmag_i_1' not in columns, see column headers: "+ str(self.columns)
		assert 'col3' in self.columns, "'col3' not in columns, see column headers: "+ str(self.columns)

		# total_bitmasks = self.data['col3']
		# if bitmask_[0] != None:
		# 	bitmask_ = bitmask_[0]

		# 	# bitmask_cut = []
		# 	# for bit_test in total_bitmasks:
		# 	# 	for mask_ in bitmask_:
		# 	# 		if (mask_ & bit_test == mask_):
		# 	# 			bitmask_cut.append(False)
		# 	# 			break
		# 	# 		if mask_ == bitmask_[-1]:
		# 	# 			bitmask_cut.append(True)

		# 	bitmask_cut = [True]*len(total_bitmasks)
		# 	for i in np.arange(0,len(bitmask_)):
		# 		# construct bitmask cut
		# 	    bitmask_cut &= np.where(bitmask_[i] & total_bitmasks == bitmask_[i], False, True)

		# 	assert len(bitmask_cut) == len(total_bitmasks), "bitmask testing broken"
		# 	bitmask_cut = np.array(bitmask_cut)
		# 	self.data = self.data[bitmask_cut]
		# 	print('bitmask cut: \t', np.unique(bitmask_cut))

		# Remove duplicates in RA/DEC:
		# coordStrings = ['RA_1', 'DEC_1']
		# for i, col in enumerate(coordStrings):
		# 	coords = self.data[col]
		# 	uniqCoords = np.unique(coords, return_inverse=True, return_counts=True)
		# 	inverse = uniqCoords[1]
		# 	count = uniqCoords[2]
		# 	orderedCount = count[inverse]
		# 	duplicateCut = orderedCount == 1
		# 	self.data = self.data[duplicateCut]
		# 	print('Removed %s duplicates in %s' % ((len(duplicateCut)-len(self.data)), col[:-2]))

		pgm = self.data['pgm']
		pgm_cut = np.array((pgm > pgm_))
		self.data = self.data[pgm_cut]
		print('pgm cut: \t', np.unique(pgm_cut))

		self.pre_count = len(self.data)
		z = self.data['z_1_1']
		self.data = self.data[(z >= 0.02)]	# define minimum redshift
		z = self.data['z_1_1']
		self.pre_z = z
		colour = self.data['absmag_g_1'] - self.data['absmag_i_1']
		total_bitmasks = self.data['col3']

		# define colour, redshift & bitmask cuts
		if colour_ != None:
			red_cut = np.array((colour > colour_)) 
			# larger (B-V) <-> 'redder' colour
			blue_cut = np.invert(red_cut)
		else:
			red_cut = np.array([True]*len(self.data))
			blue_cut = red_cut
			print('Red catalog == Blue catalog')
		print('c cut: \t', colour_, np.unique(red_cut))
		if z_ != None:
			z_cut = np.array((z > z_)) # HIGH-Z
			z_cut_r = np.invert(z_cut) # LOW-Z
		else:
			z_cut = np.array([True]*len(self.data))
			z_cut_r = z_cut
			print('highZ catalog == lowZ catalog')
		print('z cut: \t', z_, np.unique(z_cut))

		if bitmask_[0] != None:
			bitmask_ = bitmask_[0]

			# bitmask_cut = []
			# for bit_test in total_bitmasks:
			# 	for mask_ in bitmask_:
			# 		if (mask_ & bit_test == mask_):
			# 			bitmask_cut.append(False)
			# 			break
			# 		if mask_ == bitmask_[-1]:
			# 			bitmask_cut.append(True)

			bitmask_cut = [True]*len(total_bitmasks)
			for i in np.arange(0,len(bitmask_)):
				# construct bitmask cut
			    bitmask_cut &= np.where(bitmask_[i] & total_bitmasks == bitmask_[i], False, True)

			assert len(bitmask_cut) == len(total_bitmasks), "bitmask testing broken"
			bitmask_cut = np.array(bitmask_cut)
			# self.data = self.data[bitmask_cut]
			print('bitmask cut: \t', np.unique(bitmask_cut))

		# apply cuts
		self.highz_R = self.data[(z_cut & red_cut & bitmask_cut)]
		self.highz_B = self.data[(z_cut & blue_cut & bitmask_cut)]
		self.lowz_R = self.data[(z_cut_r & red_cut & bitmask_cut)]
		self.lowz_B = self.data[(z_cut_r & blue_cut & bitmask_cut)]
		self.highz = self.data[z_cut]
		self.lowz = self.data[z_cut_r]

		self.samplecounts = [len(self.highz_R), len(self.highz_B),
								len(self.lowz_R), len(self.lowz_B),
								len(self.highz), len(self.lowz)]

		[print('# objects %s: \t'%self.labels[i], v) for i, v in enumerate(self.samplecounts)]
		print('total shapes: \t%s'%np.sum(self.samplecounts[:4]))
		print('total density: \t%s'%np.sum(self.samplecounts[4:]))

		# construct sets of filenames, counts, & IDs for wcorr-calls
		self.wcorr_combos = [
		[self.labels[4]+'.asc', self.samplecounts[4], self.labels[0]+'.asc', self.samplecounts[0], 'hiZ_vs_hiZ_R'],
		[self.labels[4]+'.asc', self.samplecounts[4], self.labels[1]+'.asc', self.samplecounts[1], 'hiZ_vs_hiZ_B'],
		[self.labels[5]+'.asc', self.samplecounts[5], self.labels[2]+'.asc', self.samplecounts[2], 'loZ_vs_loZ_R'],
		[self.labels[5]+'.asc', self.samplecounts[5], self.labels[3]+'.asc', self.samplecounts[3], 'loZ_vs_loZ_B']
		]

	def cut_columns(self, subsample, h): 
		"""""
		take subsample data 
		& isolate columns for wcorr

		"""""
		assert 'RA_1' in self.columns, "'RA_1' not in columns, see column headers: "+ str(self.columns)
		assert 'DEC_1' in self.columns, "'DEC_1' not in columns, see column headers: "+ str(self.columns)
		assert 'Z_1_1' in self.columns, "'Z_1_1' not in columns, see column headers: "+ str(self.columns)
		assert 'e1c' in self.columns, "'e1c' not in columns, see column headers: "+ str(self.columns)
		assert 'e2c' in self.columns, "'e2c' not in columns, see column headers: "+ str(self.columns)

		table = subsample
		RA = np.deg2rad(table['RA_1'])
		DEC = np.deg2rad(table['DEC_1'])
		Z = table['Z_1_1']
		e1 = table['e1c']/table['pgm']
		e2 = table['e2c']/table['pgm']
		e2 *= -1 # for RA increasing leftward, c.f. x-axis increasing rightward
		e_weight = np.array([1]*len(table))

		comov = Planck13.comoving_distance(Z)
		comov *= h

		new_table = np.column_stack((RA, DEC, comov, e1, e2, e_weight))
		return new_table

	def save_tables(self, new_table, outfile_root_, label, z_cut, c_cut, notes):
		"""""
		save subsample tables to ascii

		"""""
		if notes != None:
			outfile_root_ += '_%s'%notes
		if (z_cut != None) & (c_cut != None):
			outfile_root = outfile_root_ + "_z" + str(z_cut) + "_c" + str(c_cut)
		elif z_cut != None:
			outfile_root = outfile_root_ + "_z" + str(z_cut)
		else:
			outfile_root = outfile_root_

		self.new_root = outfile_root

		if not isdir(outfile_root):
			mkdir(outfile_root)
		ascii.write(new_table, join(outfile_root, label + ".asc"), names=['#RA/rad', '#DEC/rad', '#comov_dist/Mpc/h', '#e1', '#e2', '#e_weight'])
		sample_no = str(label) + " # objects: " + str(len(new_table))
		return sample_no

	def prep_wcorr(self, files_path, wcorr_combos, rp_bins, rp_lims, los_bins, los_lim, nproc, large_pi, out_sh):

		shell_script = [
		'#!/bin/tcsh',
		'#PBS -q compute',
		'#PBS -N %s'%out_sh,
		'#PBS -l nodes=1',
		'#PBS -l walltime=120:00:00',
		'#PBS -l mem=50gb',
		'',
		'module load dev_tools/nov2014/python-anaconda',
		'',
		'cd $PBS_O_WORKDIR',
		'',
		'setenv LD_LIBRARY_PATH ${LD_LIBRARY_PATH}:/share/splinter/hj/PhD/CosmoFisherForecast/bjutils/lib/',
		'',
		'date']

		for combo in wcorr_combos: 
		# write 4x wcorr-calls to .sh script, & another 4x if largePi
			shell_script.append('')
			outfile = combo[4]
			shell_script.append('/share/splinter/hj/PhD/CosmoFisherForecast/obstools/wcorr %s %s %s %s %s %s %s %s %s %s %s %s 0 0' %	(files_path, combo[0], combo[1], combo[2], combo[3], rp_bins, rp_lims[0], rp_lims[1], los_bins, los_lim, outfile, nproc)
				)
			shell_script.append('')
			if large_pi:				
				outfile += '_largePi'
				shell_script.append('')
				shell_script.append('/share/splinter/hj/PhD/CosmoFisherForecast/obstools/wcorr %s %s %s %s %s %s %s %s %s %s %s %s %s 0' %	(files_path, combo[0], combo[1], combo[2], combo[3], rp_bins, rp_lims[0], rp_lims[1], los_bins, los_lim, outfile, nproc, large_pi)
				)
				shell_script.append('')

		shell_script.append('date')

		File = join(files_path, '%s.sh'%out_sh)
		Write = open(File, 'w')
		Text = '\n'.join(shell_script)
		Write.write(str(Text))
		Write.close()

		wcorr_spec = [] # record parameters used in wcorr
		wcorr_spec.append('Comoving transverse separation r_p: %s - %s Mpc/h in %s log-spaced bins'%(rp_lims[0], rp_lims[1], rp_bins))
		wcorr_spec.append('Comoving line-of-sight separation \Pi: %s - %s Mpc/h in %s bins'%(los_lim*(-1), los_lim, los_bins))
		wcorr_spec.append('No. processors: %s'%nproc)
		wcorr_spec.append('Large-Pi systematics testing: %s'%large_pi)

		File = join(files_path, 'wcorr_spec')
		Write = open(File, 'w')
		Text = '\n'.join(wcorr_spec)
		Write.write(str(Text))
		Write.close()

	def plot_wcorr(self, files_path, wcorrIDs):
		# 0 = hiZ_Red
		# 1 = hiZ_Blue
		# 2 = loZ_Red
		# 3 = loZ_Blue
		wcorrOutputs = []
		rand_wcorrOutputs = []
		for item in wcorrIDs:
			# construct filenames.dat
			wcorrOutputs.append('%s'%join(files_path, ('wcorr_' + item + '.dat')))
			rand_wcorrOutputs.append('%s'%join(files_path, ('wcorr_rand_' + item + '.dat')))
		realData = []
		randData = []
		wgplus = []
		wgcross = []
		wgerr = []
		rand_wgplus = []
		rand_wgcross = []
		rand_wgerr = []
		for i, path in enumerate(wcorrOutputs):
			realData.append(np.loadtxt(path))
			randData.append(np.loadtxt(rand_wcorrOutputs[i]))
			# subtract randoms from +/x
			realData[i][:,3] -= randData[i][:,3]
			realData[i][:,4] -= randData[i][:,4]
			realErr = realData[i][:,6]
			randErr = randData[i][:,6]
			# propagate errors
			propgErrs = np.sqrt((realErr**2) + (randErr**2))
			wgplus.append(realData[i][:,3])
			wgcross.append(realData[i][:,4])
			wgerr.append(propgErrs)
			rand_wgplus.append(randData[i][:,3])
			rand_wgcross.append(randData[i][:,4])
			rand_wgerr.append(randData[i][:,6])
		r_p = realData[0][:,0]
		x = np.linspace(0, r_p.max()*1.8)
		# plt.ioff()
		dataPoints = [[wgplus, wgcross, wgerr], [rand_wgplus, rand_wgcross, rand_wgerr]]
		prefix = ['','rand_']
		for j, set_ in enumerate(dataPoints): 
		# plot/save random-subtracted-reals, AND randoms
			f, axarr = plt.subplots(2, 2, sharex=True, sharey=True, figsize=(15,10))
			f.subplots_adjust(hspace=0, wspace=0)
			axarr[0,0].errorbar(r_p, set_[0][0], yerr=set_[2][0],
								elinewidth=2, color='r', capsize=0,
								label='w(g+)')
			axarr[0,0].errorbar(r_p, set_[1][0], yerr=set_[2][0],
								elinewidth=2, color='g', capsize=0,
								label='w(gx)', alpha=0.5)
			axarr[0,1].errorbar(r_p, set_[0][1], yerr=set_[2][1],
								elinewidth=2, color='b', capsize=0,
								label='w(g+)')
			axarr[0,1].errorbar(r_p, set_[1][1], yerr=set_[2][1],
								elinewidth=2, color='g', capsize=0,
								label='w(gx)', alpha=0.5)
			axarr[1,0].errorbar(r_p, set_[0][2], yerr=set_[2][2],
								elinewidth=2, color='r', capsize=0,
								label='w(g+)')
			axarr[1,0].errorbar(r_p, set_[1][2], yerr=set_[2][2],
								elinewidth=2, color='g', capsize=0,
								label='w(gx)', alpha=0.5)
			axarr[1,1].errorbar(r_p, set_[0][3], yerr=set_[2][3],
								elinewidth=2, color='b', capsize=0,
								label='w(g+)')
			axarr[1,1].errorbar(r_p, set_[1][3], yerr=set_[2][3],
								elinewidth=2, color='g', capsize=0,
								label='w(gx)', alpha=0.5)
			arr_ind = [(0,0), (0,1), (1,0), (1,1)]
			for i, ind in enumerate(arr_ind):
				a = axarr[ind]
				a.set_xscale('log')
				a.set_xlim(0.25,70)
				a.set_ylim(-0.5,0.4)
				a.plot(x, [0]*len(x), lw=2, ls='--', color='c')
				# a.set_xlabel('Comoving transverse separation (Mpc/h)')
				# a.set_ylabel('Correlations')
				a.set_title('%s'%wcorrIDs[i], fontsize=12)
				a.legend(loc='upper right')
				a.grid()
			plt.setp([a.get_xticklabels() for a in axarr[0, :]], visible=False)
			plt.setp([a.get_yticklabels() for a in axarr[:, 1]], visible=False)
			# plt.setp([a.get_xlabels() for a in axarr[:, 1]], visible=False)
			axarr[1,0].set_xlabel('Comoving transverse separation (Mpc/h)')
			axarr[1,0].set_ylabel('Correlations')

			plotsDir = join(files_path, 'Plots')
			if not isdir(plotsDir):
				mkdir(plotsDir)

			if 'largePi' in wcorrOutputs[0]:
				outimg = join(plotsDir, '%swcorr_plots_largePi.pdf'%prefix[j])
			else:
				outimg = join(plotsDir, '%swcorr_plots.pdf'%prefix[j])
			f.savefig(outimg)

		return wcorrOutputs

	def chiFunc(self, y):
		return chi2.pdf(y, 10)

	def normFunc(self, y):
		return stat.norm(0,1).pdf(y)

	def chi2(self, path2data, expec):
		filesList = np.array(listdir(path2data))
		datCut = np.array([i.endswith('.dat') for i in filesList])
		dataList = filesList[datCut]
		dataArr = [np.loadtxt(join(path2data, i)) for i in dataList]
		dataArr = np.array([[i[:,3], i[:,4], i[:,6]] for i in dataArr])
		randCut = np.array([i.startswith('wcorr_rand') for i in dataList])
		realCut = np.invert(randCut)
		randData = dataArr[randCut]
		dataArr = dataArr[realCut]
		pVals = []
		chiSqs = []
		xSigma = []

		for j, data in enumerate(dataArr):
			plus = data[0]-randData[j][0]
			cross = data[1]-randData[j][1]
			err = np.sqrt(data[2]**2+randData[j][2]**2)#[1]*df
			plusChi_i = [((v-expec)/err[i])**2 for i, v in enumerate(plus)]
			crossChi_i = [((v-expec)/err[i])**2 for i, v in enumerate(cross)]
			chiSq_pl = np.sum(plusChi_i)
			chiSq_cr = np.sum(crossChi_i)
			intgrl_pl = scint.quad(self.chiFunc, chiSq_pl, np.inf)
			intgrl_cr = scint.quad(self.chiFunc, chiSq_cr, np.inf)
			pVal_pl = intgrl_pl[0]
			pVal_cr = intgrl_cr[0]
			pVal = ['%.5f'%pVal_pl, '%.5f'%pVal_cr]
			chiSq = ['%.5f'%chiSq_pl, '%.5f'%chiSq_cr]
			chiSqs.append(chiSq)
			pVals.append(pVal)

			xSigs = []
			for p in [pVal_pl, pVal_cr]:
				x = p/10
				int_x = scint.quad(self.normFunc,x,np.inf)
				gauss_p = 1-(2*int_x[0])
				gaussOver_p = gauss_p/p
				while abs(1-gaussOver_p) > 0.01:
					x *= 1.01
					int_x = scint.quad(self.normFunc,x,np.inf)
					gauss_p = 1-(2*int_x[0])
					gaussOver_p = gauss_p/p
				else:
					xSigs.append(['%.5f'%x, '%.5f'%gauss_p])
			xSigma.append(xSigs)

		for l in [pVals,chiSqs,xSigma]:
			l = np.array(l)
		chi2Stats = zip(dataList,chiSq[:,0],pVal[:,0],xSigma[:,0,0],xSigma[:,0,1],chiSq[:,1],pVal[:,1],xSigma[:,1,0],xSigma[:,1,1])
		fl = open(join(path2data, 'chisquares.csv'))
		writer = csv.writer(fl)
		writer.writerow(['data','chi^2(plus)','p-val','x-sigma','(gauss-p)','chi^2(cross)','p-val','x-sigma','(gauss-p)'])
		for vals in chi2Stats:
			writer.writerow(vals)
		fl.close()

		return None


class RandomCatalogue(RealCatalogue):

	def __init__(self, path):
		"""""
		read-in catalogue

		"""""
		self.path = path
		hdulist = fits.open(path)
		self.data = hdulist[1].data
		self.columns = hdulist[1].columns.names
		self.labels = ['rand_highZ', 'rand_lowZ']

	def cut_data(self, z_, len_reals, z_reals): 
		"""""
		cut catalogue into redshift subsamples
		
		"""""
		assert 'RA' in self.columns, "'RA' not in columns, see column headers: "+ str(self.columns)
		assert 'DEC' in self.columns, "'DEC' not in columns, see column headers: "+ str(self.columns)
		assert 'Z' in self.columns, "'Z' not in columns, see column headers: "+ str(self.columns)
		assert 'RAND_NUM' in self.columns, "'RAND_NUM' not in columns, see column headers: "+ str(self.columns)

		# cut down (randomly) for correlation
		randNum = self.data['RAND_NUM']
		current = len(randNum)
		target = len_reals*10
		fraction = target/current
		randCut = (randNum > fraction) & (randNum <= 2*fraction)
		self.data = self.data[randCut]
		print('Randoms cut down from %s objects to %s' % (len(randNum), len(self.data)))

		# Remove duplicates in RA/DEC:
		# coordStrings = ['RA', 'DEC']
		# for i, col in enumerate(coordStrings):
		# 	coords = self.data[col]
		# 	uniqCoords = np.unique(coords, return_inverse=True, return_counts=True)
		# 	inverse = uniqCoords[1]
		# 	count = uniqCoords[2]
		# 	orderedCount = count[inverse]
		# 	duplicateCut = orderedCount == 1
		# 	self.data = self.data[duplicateCut]
		# 	print('Removed %s duplicates in %s' % ((len(duplicateCut)-len(self.data)), col))

		z = self.data['z']
		pre_z_cut = (z >= z_reals.min()) & (z <= z_reals.max())
		self.data = self.data[pre_z_cut]
		z = self.data['z']

		if z_ != None:
			z_cut = np.array((z > z_))
			z_cut_r = np.invert(z_cut)
		else:
			z_cut = np.array([True]*len(self.data))
			z_cut_r = z_cut
			print('highZ catalog == lowZ catalog')

		self.highz = self.data[z_cut]
		self.lowz = self.data[z_cut_r]
		self.samplecounts = [len(self.highz), len(self.lowz)]

		[print('# objects %s: \t'%self.labels[i], v) for i, v in enumerate(self.samplecounts)]

	def cut_columns(self, subsample, h): 
		"""""
		take subsample data 
		& isolate columns for wcorr

		"""""
		assert 'RA' in self.columns, "'RA' not in columns, see column headers: "+ str(self.columns)
		assert 'DEC' in self.columns, "'DEC' not in columns, see column headers: "+ str(self.columns)
		assert 'Z' in self.columns, "'Z' not in columns, see column headers: "+ str(self.columns)

		table = subsample
		RA = np.deg2rad(table['RA'])
		DEC = np.deg2rad(table['DEC'])
		Z = table['Z']
		e1 = 2*np.random.random(len(table))
		e2 = e1
		e2 *= -1 # for RA increasing leftward, c.f. x-axis increasing rightward
		e_weight = np.array([1]*len(table))

		comov = Planck13.comoving_distance(Z)
		comov *= h

		new_table = np.column_stack((RA, DEC, comov, e1, e2, e_weight))
		return new_table

if __name__ == "__main__":
	parser = argparse.ArgumentParser(fromfile_prefix_chars='@')
	parser.add_argument(
		'Catalog',
		help='full path of REAL catalogue to be sampled into ascii table(s)')
	parser.add_argument(
		'Path',
		help='full path of destination directory for subsample ascii catalogues, where further directories will be created and appended with "_z_<z_cut>_c_<colour_cut>" if applicable. ***If using -plotNow, give full path to wcorr output directory***')
	parser.add_argument(
		'-Random',
		help="optional; path to RANDOM catalogue to be correspondingly sampled",
		default=None)
	parser.add_argument(
		'-pgm_cut',
		type=np.float32,
		help='shear polarizability cut, defaults to 0.1',
		default=0.1)
	parser.add_argument(
		'-zCut',
		type=np.float32,
		help='lowZ vs. highZ redshift threshold, between 0 - 0.5. Omit for no redshift cut',
		default=None)
	parser.add_argument(
		'-cCut',
		type=np.float32,
		help='red vs. blue colour threshold, between 0 - 1.4 (meaningfully, between approx 0.7 - 1.1). Omit for no colour cut',
		default=None)
	parser.add_argument(
		'-bitmaskCut',
		nargs='*',
		type=int,
		help='list of bitmask IDs (powers of 2) to exclude from catalogue, eg. text-file w/ ROWS; 4, 16, 4096...etc...',
		default=None)
	parser.add_argument(
		'-H',
		type=np.float32,
		help='reduced Planck constant, defaults to 0.7',
		default=0.7)
	parser.add_argument(
		'-rpBins',
		type=int,
		help='specify no. of (log-spaced) bins in comoving transverse separation r_p (Mpc/h), for measurement of density-shape correlations. Defaults to 10',
		default=10)
	parser.add_argument(
		'-rpLims',
		nargs=2,
		type=np.float32,
		help='specify upper & lower (2 args, space-separated) limit in comoving transverse separation r_p (Mpc/h), for measurement of density-shape correlations. Defaults to 0.3, 60 Mpc/h',
		default=[0.3, 60])
	parser.add_argument(
		'-losBins',
		type=int,
		help='specify no. of bins (positive branch) in comoving line-of-sight separation \Pi (Mpc/h), for measurement of density-shape correlations. Defaults to 15; i.e. 30 bins in total',
		default=15)
	parser.add_argument(
		'-losLim',
		type=np.float32,
		help='specify cut-off in comoving line-of-sight separation \Pi (Mpc/h), for measurement of density-shape correlations. Defaults to 60; s/t range is -60 to 60 Mpc/h',
		default=60)
	parser.add_argument(
		'-nproc',
		type=int,
		help='no. processors to be used in correlation measurement, default = 12',
		default=12)
	parser.add_argument(
		'-largePi',
		type=int,
		choices=[0,1],
		help='specify regular (0) or regular + large-Pi systematics tests (1), defaults to 1',
		default=0)
	parser.add_argument(
		'-wcorr',
		type=int,
		choices=[0,1],
		help='initiate wcorr density-shape correlation measurements (1), or not (0), defaults to 0',
		default=0)
	parser.add_argument(
		'-notes',
		help='notes on any changed wcorr parameters, for appendage to directory name',
		default=None)
	parser.add_argument(
		'-plot',
		help='plot data & save figure after correlations completed (1), or not (0), defaults to 1',
		default=1)
	parser.add_argument(
		'-plotNow',
		help='plot ALREADY EXISTING correlation data (1), having given arg="Path" as the path to the .dat files (Catalog arg must still be path of readable .fits catalog). Bypasses all other sampling functions. Defaults to 0',
		default=0)
	parser.add_argument(
		'-chiSqu',
		help='calc chi^2 stats for ALREADY EXISTING correlation data (1), having given arg="Path" as the path to the .dat files (Catalog arg must still be path of readable .fits catalog). Bypasses all other sampling functions. Defaults to 0',
		default=0)
	parser.add_argument(
		'-expec',
		help='expectation values for chi^2 statistics, defaults to zeros at all points',
		default=0)
	args = parser.parse_args()

	catalog = RealCatalogue(args.Catalog)

	if args.plotNow:
		# plot .dat files, returning filename-list
		wcorrOuts = catalog.plot_wcorr(args.Path, catalog.wcorrLabels)
		largePi_outs = [basename(normpath(out[:-4] + '_largePi.dat')) for out in wcorrOuts]
		# check for largePi .dat files
		isIn = [i in listdir(args.Path) for i in largePi_outs]
		uniq = np.unique(isIn)
		if uniq.all() == True:
			IDs = [outs[6:-4] for outs in largePi_outs]
			a = catalog.plot_wcorr(args.Path, IDs)

	# 	if args.chiSqu:
	# 		# calculate chi^2 statistics & save to csv
	# 		catalog.chi2(args.Path, args.expec)
	# 	sys.exit()

	# if args.chiSqu:
	# 	# calculate chi^2 statistics & save to csv
	# 	catalog.chi2(args.Path, args.expec)
	# 	sys.exit()

	catalog.cut_data(args.pgm_cut, args.zCut, args.cCut, args.bitmaskCut)
	samples = [catalog.highz_R, catalog.highz_B, 									catalog.lowz_R, catalog.lowz_B,										catalog.highz, catalog.lowz]
	cuts = 'z-cut: ' + str(args.zCut) + ', colour-cut (g-i): ' + str(args.cCut)
	sample_numbers = [cuts]
	outfile_root = join(args.Path, 'Wcorr')

	for i, sample in enumerate(samples):
		new_table = catalog.cut_columns(sample, args.H)
		sample_num = catalog.save_tables(new_table, outfile_root, catalog.labels[i], args.zCut, args.cCut, args.notes)
		sample_numbers.append(sample_num)

	File = join(catalog.new_root, 'Sample_popns')
	Write = open(File, "w")
	Text = "\n".join(sample_numbers)
	Write.write(str(Text))
	Write.close()

	catalog.prep_wcorr(catalog.new_root, catalog.wcorr_combos, args.rpBins, args.rpLims, args.losBins, args.losLim, args.nproc, args.largePi, 'real_wcorr')

	if args.wcorr:
		list_dir = np.array(listdir(catalog.new_root))
		shells = np.array([i.endswith('.sh') for i in list_dir])
		list_dir = list_dir[shells]
		[os.system('qsub '+ join(catalog.new_root, shell)) for shell in list_dir]

	if args.Random != None:
		catalog2 = RandomCatalogue(args.Random)
		catalog2.cut_data(args.zCut, catalog.pre_count, catalog.pre_z)
		samples = [catalog2.highz, catalog2.lowz]
		cuts = 'z-cut: ' + str(args.zCut)
		sample_numbers = [cuts]

		for i, sample in enumerate(samples):
			new_table = catalog2.cut_columns(sample, args.H)
			sample_num = catalog2.save_tables(new_table, outfile_root, catalog2.labels[i], args.zCut, args.cCut, args.notes)
			sample_numbers.append(sample_num)

		File = join(catalog.new_root, 'Sample_rand_popns')
		Write = open(File, "w")
		Text = "\n".join(sample_numbers)
		Write.write(str(Text))
		Write.close()

		rand_combos = [
		[catalog2.labels[0]+'.asc', catalog2.samplecounts[0], catalog.labels[0]+'.asc', catalog.samplecounts[0], 'rand_hiZ_vs_hiZ_R'],
		[catalog2.labels[0]+'.asc', catalog2.samplecounts[0], catalog.labels[1]+'.asc', catalog.samplecounts[1], 'rand_hiZ_vs_hiZ_B'],
		[catalog2.labels[1]+'.asc', catalog2.samplecounts[1], catalog.labels[2]+'.asc', catalog.samplecounts[2], 'rand_loZ_vs_loZ_R'],
		[catalog2.labels[1]+'.asc', catalog2.samplecounts[1], catalog.labels[3]+'.asc', catalog.samplecounts[3], 'rand_loZ_vs_loZ_B']
		]

		catalog2.prep_wcorr(catalog.new_root, rand_combos, args.rpBins, args.rpLims, args.losBins, args.losLim, args.nproc, args.largePi, 'rand_wcorr')

		if args.plot:
			with open(join(catalog.new_root, 'rand_wcorr.sh'), 'a') as script:
				script.write(
					'\npython /share/splinter/hj/PhD/catalog_sampler.py %s %s -plotNow 1 -chiSqu 1'%(args.Catalog, catalog.new_root)
					)
				script.write('\n')

		if args.wcorr:
			list_dir = np.array(listdir(catalog.new_root))
			shells = np.array([i.endswith('.sh') for i in list_dir])
			r_shells = np.array([i.startswith('rand') for i in list_dir])
			list_dir = list_dir[(shells&r_shells)]
			[os.system('qsub '+ join(catalog.new_root, shell)) for shell in list_dir]

























