# coding: utf-8
from __future__ import print_function, division ; import numpy as np ; import matplotlib.pylab as plt ; import matplotlib.lines as ml
import os
from os import listdir, mkdir
from os.path import join, isdir, basename, normpath, dirname
import sys
sys.path.append("/share/splinter/hj/PhD/")
from jackknife3d import betwixt

def read_cats(sample_path, jk_number=0, shapes_patch=0):
	dens_path = basename(normpath(sample_path)) + '_UnMasked'
	if shapes_patch:
		shapes = np.loadtxt( join( sample_path, sample_path + '%spatch.asc'%str(jk_number).zfill(3) ))
	else:
		shapes = np.loadtxt( join( sample_path, 'JKsamples', 'JKsample%s.asc'%str(jk_number).zfill(3) ))
	densities = np.loadtxt( join( dens_path, 'JKsamples', 'JKsample%s.asc'%str(jk_number).zfill(3) ))
	randoms = np.loadtxt( join( dens_path, 'JKsamples', 'rand_JKsample%s.asc'%str(jk_number).zfill(3) ))
	return shapes, densities, randoms

def plot3dhist(shapes, densities=None, randoms=None, sdss=0, legend=1):
    f, ax = plt.subplots(3, figsize=(10,8)) ; plt.subplots_adjust(hspace=0.35)
    for x in range(3):
        dim = ('ra','dec','z')[x]
	ax[x].set_xlabel(dim)
        bins = (np.linspace(shapes.T[0].min(),shapes.T[0].max(),101), np.linspace(shapes.T[1].min(),shapes.T[1].max(),101), np.linspace(shapes.T[2].min(),shapes.T[2].max(),101))[x]
        h_ra, e_ra, p = ax[x].hist(shapes.T[x], bins=bins, stacked=True, normed=True, color='r', histtype='step', lw=2, alpha=0.7, label='shapes %s'%dim)
	try:
		h_dec, e_dec, p = ax[x].hist(densities.T[x], bins=bins, stacked=True, normed=True, color='b', histtype='step', lw=2, alpha=0.7, label='density %s'%dim)
		h_z, e_z, p = ax[x].hist(randoms.T[x], bins=bins, stacked=True, normed=True, color='g', histtype='step', lw=2, alpha=0.7, label='random %s'%dim)
		ax[x].set_ylim( [h_ra, h_dec, h_z][x].min(), [h_ra, h_dec, h_z][x].max()*1.1)
	except AttributeError:
		legend = 0
		continue
    if legend:
	h = [ml.Line2D([], [], ls='-', lw=2, c=('r','b','g')[i], label=('shapes','densities','randoms')[i]) for i in range(3)]
	if sdss:
	    ax[0].legend(handles=h, loc='upper left', title='normalised hists', fontsize=8, ncol=2)
	else:
	    ax[-1].legend(handles=h, loc='best', title='normalised hists', fontsize=8, ncol=2)
	
    
def plot3dims(shapes, densities=None, randoms=None, coordinate_bounds=None):
    f, ax = plt.subplots(3, figsize=(12,10)) ; plt.subplots_adjust(hspace=0.35)
    if coordinate_bounds!=None:
	xmin, xmax, ymin, ymax, zmin, zmax = coordinate_bounds
	boxcutter = betwixt( coordinate_bounds )
	print("x-range: %.2f - %.2f \n" % (xmin, xmax),
		"y-range: %.2f - %.2f \n" % (ymin, ymax),
		"z-range: %.2f - %.2f \n" % (zmin, zmax))
	print('shapes count: ', sum( boxcutter( shapes.T[0], shapes.T[1], shapes.T[2] ) ) )		## <-- UGLY :( ## re-write betwixt() output (& all fns that call it) to take an array
	shapes = shapes[ boxcutter( shapes.T[0], shapes.T[1], shapes.T[2] ) ]
	try:
		print('density count: ', sum( boxcutter( densities.T[0], densities.T[1], densities.T[2] ) ) ) 
		print('randoms count: ', sum( boxcutter( randoms.T[0], randoms.T[1], randoms.T[2] ) ) )
		densities = densities[ boxcutter( densities.T[0], densities.T[1], densities.T[2] ) ]
		randoms = randoms[ boxcutter( randoms.T[0], randoms.T[1], randoms.T[2] ) ]
	except AttributeError:
		print('not plotting shapes vs. densities vs. randoms')

    for i, x in enumerate([(0, 1), (1, 2), (2, 0)]):
	try:
		ax[i].plot(randoms.T[x[0]], randoms.T[x[1]], c='g', marker='.', markeredgecolor='g', ms=2, alpha=0.3, ls='')
		ax[i].plot(densities.T[x[0]], densities.T[x[1]], c='b', marker='.', markeredgecolor='b', ms=2, alpha=0.4, ls='')
	except AttributeError:
		print('not plotting shapes vs. densities vs. randoms')
        ax[i].plot(shapes.T[x[0]], shapes.T[x[1]], c='r', marker='.', markeredgecolor='r', ms=2, alpha=0.4, ls='')
        ax[i].set_xlabel(('ra','dec','$\chi$')[x[0]])
        ax[i].set_ylabel(('ra','dec','$\chi$')[x[1]])
        ax[i].grid(which='both', ls='-', alpha=0.2)
    #h = [ml.Line2D([], [], ls='-', lw=2, c=('r','b','g')[i], label=('shapes','densities','randoms')[i]) for i in range(3)]
    
