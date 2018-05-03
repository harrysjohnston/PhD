# coding: utf-8
from __future__ import print_function, division ; import numpy as np ; import matplotlib.pylab as plt ; import matplotlib.lines as ml
import os
from os import listdir, mkdir
from os.path import join, isdir, basename, normpath, dirname
import sys
sys.path.append("/share/splinter/hj/PhD/")
from jackknife3d import betwixt

def read_ncats(*paths):

	catalog_array = []
	for path in paths:
		print(path)
		cat = np.loadtxt(path)
		catalog_array.append(cat)
	#print('SDSS ?')

	return dict(zip(paths, np.array(catalog_array)))

#def read_cats(sample_path, jk_number=0, shapes_patch=0):
#	dens_path = basename(normpath(sample_path)) + '_UnMasked'
#	if shapes_patch:
#		shapes = np.loadtxt( join( sample_path, sample_path + '%spatch.asc'%str(jk_number).zfill(3) ))
#	else:
#		shapes = np.loadtxt( join( sample_path, 'JKsamples', 'JKsample%s.asc'%str(jk_number).zfill(3) ))
#	densities = np.loadtxt( join( dens_path, 'JKsamples', 'JKsample%s.asc'%str(jk_number).zfill(3) ))
#	randoms = np.loadtxt( join( dens_path, 'JKsamples', 'rand_JKsample%s.asc'%str(jk_number).zfill(3) ))
#	return shapes, densities, randoms

def plot3dhist(sdss=0, legend=1, binning=100, **cats):
    keys = cats.keys()
    keys.sort()
    keys = keys[::-1]
    f, ax = plt.subplots(3, figsize=(10,8)) ; plt.subplots_adjust(hspace=0.35)
    cspace = np.linspace(0.1,0.9, len(keys))
    cmap = plt.cm.brg(cspace)
    for x in range(3):
        dim = ('ra','dec','z')[x]
	ax[x].set_xlabel(dim)
	for i, k in enumerate(keys):
        	bins = ( np.linspace(cats[k].T[0].min(), cats[k].T[0].max(), binning + 1), np.linspace(cats[k].T[1].min(), cats[k].T[1].max(), binning+ 1), np.linspace(cats[k].T[2].min(), cats[k].T[2].max(), binning + 1) ) [x]
		lab = basename( normpath(k) )
		h_x, e_x, p = ax[x].hist(cats[k].T[x], bins=bins, stacked=True, normed=True, color=cmap[i], histtype='step', lw=2, alpha=0.7, label=lab)
		ax[x].set_ylim( h_x.min(), h_x.max()*1.1 )
	ax[x].grid(which='both', ls='-', alpha=0.2)

    if legend:
	h = [ ml.Line2D( [], [], ls='-', lw=2, c=cmap[i], label=basename( normpath(k) ) ) for i, k in  enumerate(keys) ]
	if sdss:
	    ax[1].legend(handles=h, loc='upper center', fontsize=8, ncol=len(keys)//2)
	else:
	    ax[1].legend(handles=h, loc='lower center', fontsize=8, ncol=len(keys)//2)
	plt.tight_layout()
	return ax
	
    
def plot3dims(sdss=0, legend=1, **cats):
    keys = cats.keys()
    keys.sort()
    keys = keys[::-1]
    f, ax = plt.subplots(3, figsize=(12,10)) ; plt.subplots_adjust(hspace=0.35)
    cspace = np.linspace(0.1,0.9, len(keys))
    cmap = plt.cm.brg(cspace)
#    if coordinate_bounds!=None:
#	xmin, xmax, ymin, ymax, zmin, zmax = coordinate_bounds
#	boxcutter = betwixt( coordinate_bounds )
#	print("x-range: %.2f - %.2f \n" % (xmin, xmax),
#		"y-range: %.2f - %.2f \n" % (ymin, ymax),
#		"z-range: %.2f - %.2f \n" % (zmin, zmax))
#	print('shapes count: ', sum( boxcutter( shapes.T[0], shapes.T[1], shapes.T[2] ) ) )		## <-- UGLY :( ## re-write betwixt() output (& all fns that call it) to take an array
#	shapes = shapes[ boxcutter( shapes.T[0], shapes.T[1], shapes.T[2] ) ]
#	try:
#		print('density count: ', sum( boxcutter( densities.T[0], densities.T[1], densities.T[2] ) ) ) 
#		print('randoms count: ', sum( boxcutter( randoms.T[0], randoms.T[1], randoms.T[2] ) ) )
#		densities = densities[ boxcutter( densities.T[0], densities.T[1], densities.T[2] ) ]
#		randoms = randoms[ boxcutter( randoms.T[0], randoms.T[1], randoms.T[2] ) ]
#	except AttributeError:
#		print('not plotting shapes vs. densities vs. randoms')

    for i, x in enumerate([(0, 1), (1, 2), (2, 0)]):
	for j, k in enumerate(keys):
		ax[i].plot(cats[k].T[x[0]], cats[k].T[x[1]], c=cmap[j], marker='.', ms=2, alpha=0.3, ls='')#, markeredgecolor=cmap[j])
		if (i == 0) & (j%2==0): ax[i].annotate(j, xy=(cats[k].T[x[0]].mean(), cats[k].T[x[1]].mean()), xycoords='data', fontsize=11)#, markeredgecolor=cmap[j])
        
	ax[i].set_xlabel(('ra','dec','$\chi$')[x[0]])
        ax[i].set_ylabel(('ra','dec','$\chi$')[x[1]])
        ax[i].grid(which='both', ls='-', alpha=0.2)

    if legend:
	h = [ ml.Line2D( [], [], ls='-', lw=2, c=cmap[i], label=basename( normpath(k) ) ) for i, k in  enumerate(keys) ]
	if sdss:
	    ax[0].legend(handles=h, loc='lower right', fontsize=8, ncol=len(keys)//2)
	else:
	    ax[0].legend(handles=h, loc='lower left', fontsize=8, ncol=len(keys)//2)
	plt.tight_layout()
	return ax

    
