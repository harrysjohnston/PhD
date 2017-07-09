# coding: utf-8
from __future__ import print_function, division ; from astropy.io import fits, ascii ; import numpy as np ; import matplotlib.pylab as plt
from astropy.cosmology import Planck13
import healpy as hp
import gc

def betwixt((re1,re2,de1,de2,ze1,ze2)):
	def make_cut(ra,dec,z):
		return (ra>=re1)&(ra<=re2)&(dec>=de1)&(dec<=de2)&(z>=ze1)&(z<=ze2)
	return make_cut

def resample_data(fitsdata, sample_cuts, patchside=6, do_sdss=0, do_3d=1, cube_zdepth=0.06, largePi=0, bitmaskCut=None):
	resampler = resampleTools(fitsdata, patchside, do_sdss, do_3d, sample_cuts, cube_zdepth=cube_zdepth, largePi=largePi)

	# identify masked pixel coordinates
	lostpix_coords = resampler.find_lostpixels(bitmaskCut)
	# lostpix_coords = (np.array([1e3]),np.array([1e3]))

	# define patches & their degree of masking
	patch_cuts, patch_weights, patch_idx, edges = resampler.define_edgecuts(lostpix_coords)

	# apply patch cuts to data
	patches = resampler.make_patches(fitsdata, patch_cuts)

	# filter patches straddling survey-edges
	patch_cuts, patch_weights, error_scaling = resampler.filter_patches(patches, patch_cuts, patch_weights, patch_idx, edges,lilbin_division=4,occupn_threshold=0.9)

	# create function to trim randoms down to real JK footprint
	random_cutter = resampler.find_patch_edges(patch_cuts)

	# apply sample cuts to patch-cuts
	sample_patch_cuts = []
	for s_cut in resampler.sample_cuts:
		spc = resampler.make_sample_patchcuts(patch_cuts, s_cut)
		sample_patch_cuts.append(spc)
	sample_patch_cuts = np.array(sample_patch_cuts)

	# construct array of patches, per sample
	patchData = []
	for spc in sample_patch_cuts:
		sample_patches = resampler.make_patches(fitsdata, spc)
		patchData.append(sample_patches)

	if do_3d:
		# slice patches in z, creating jackknife cubes
		cubeData = []
		for spatches in patchData:
			cubes, cube_weights = resampler.make_cubes(spatches, patch_weights, edges[2])
			cubeData.append(cubes)
	else:
		patchData = np.array(patchData)
		cubeData = patchData.copy()
		cube_weights = patch_weights.copy()
		del patchData, patch_weights

	cubeData = np.array(cubeData)

	print('\nN cubes: ', len(cubeData[0]),
			'\nmin | max cube weights: %.4f | %.4f'%(cube_weights.min(), cube_weights.max()))

	return cubeData, cube_weights, error_scaling, random_cutter

class resampleTools:
	def __init__(self, fitsdata, patchside, do_sdss, do_3d, sample_cuts, cube_zdepth=0.06, largePi=0):
		self.cols = [ ['RA_GAMA', 'DEC_GAMA', 'Z_TONRY'], ['ra', 'dec', 'z'] ] [do_sdss]
		self.ranges = [ [ (i, j, -3., 3., 0., 0.5) for i,j in [ (129.,141.), (174.,186.), (211.5,223.5) ] ], None ] [do_sdss]
		self.data = fitsdata
		self.do_sdss = do_sdss
		self.do_3d = do_3d
		self.ra_side = self.dec_side = patchside
		self.z_side = cube_zdepth
		self.sample_cuts = sample_cuts
		self.largePi = largePi

	def find_lostpixels(self, *bitmask_):
		# find coordinates of pixels lost to masking
		print('finding pixels lost to mask..')
		if not self.do_sdss:
			nside = 2048
			fullSky = 41252.96 # square degrees
			npix = hp.nside2npix(nside)
			ra = self.data['RA_GAMA']
			dec = self.data['DEC_GAMA']
			theta = np.deg2rad(90.-dec)
			phi = np.deg2rad(ra)
			pixIDs = hp.ang2pix(nside,theta,phi,nest=False)
			GKmap = np.bincount(pixIDs,minlength=npix)
			GKpix = np.where(GKmap!=0,1,0)
			# kidsBitmap = hp.read_map('/share/splinter/hj/PhD/KiDS_counts_N2048.fits', dtype=int)
			kidsBitmap = hp.read_map('./KiDS_counts_N2048.fits', dtype=int)
			bitmask_cut = [True]*len(kidsBitmap)

			print("PATCHES: CUTTING MASK!=0")
			bitmask_cut = np.where(kidsBitmap==0,True,False)

			if bitmask_[0] != None:
				bitmask_ = bitmask_[0]
				for i in range(len(bitmask_)):
					# construct specific bitmask cut
					bitmask_cut &= np.where(bitmask_[i] & kidsBitmap == bitmask_[i], False, True)

			lostpixIDs = [j for j,b in enumerate(bitmask_cut) if b==False]
			lostfromGK = np.where(GKpix[lostpixIDs]!=0,1,0) #1=pixel-lost,0=not-lost
			print('Lost npix, fraction of area: %s, %.3f'%(sum(lostfromGK),sum(lostfromGK)/sum(GKpix)))

			if lostpixIDs != []:
				thetaPhis = hp.pix2ang(nside,lostpixIDs)
				# lost pixel coords;
				lostpixra,lostpixdec = np.rad2deg(thetaPhis[1]),(90.-np.rad2deg(thetaPhis[0]))
			else:
				lostpixra,lostpixdec = (np.array([1e3]),np.array([1e3]))
			del kidsBitmap
			gc.collect()
		else:
			lostpixra,lostpixdec = (np.array([1e3]),np.array([1e3]))
		# lostpixra,lostpixdec = (np.array([1e3]),np.array([1e3]))

		lostpix_coords = (lostpixra, lostpixdec)
		return lostpix_coords

	def define_edgecuts(self, lostpix_coords):
		# given desired patch/box-sizes, divide sky into patches
		# return sets of patch-cuts for application to catalogs, with weights due to lost pixels

		gal_coords = np.column_stack((self.data[self.cols[0]], self.data[self.cols[1]], self.data[self.cols[2]]))
		ra, dec, z = gal_coords.T
		lpix_ra, lpix_dec = lostpix_coords

		print('defining cubes..')
		ra_num, dec_num, z_num = 360//self.ra_side +1, 180//self.dec_side +1, 0.5//self.z_side +1
		if self.largePi:
			z_num = 0.5//0.1 +1
		redg,dedg,zedg = np.linspace(0,360,ra_num), np.linspace(-90,90,dec_num), np.linspace(0,0.5,z_num)
		if self.ranges!=None:
			print('GAMA equatorial: jackknifing specified ranges...')
			redg = []
			dec_num = 6.//self.dec_side +1
			dedg = np.linspace(-3., 3., dec_num)
			zedg = list(np.linspace(0., 0.5, z_num))
			for edges in self.ranges:
				ra_num = (edges[1] - edges[0])//self.ra_side + 1
				redg += list(np.linspace(edges[0], edges[1], ra_num))

		redg,dedg,zedg = map(lambda x: np.array(x), [redg,dedg,zedg])
		radiff, decdiff, zdiff = (np.diff(i) for i in [redg,dedg,zedg])
		print('ra | dec | z sides: %.2f | %.2f | %.3f'%(radiff.min(),decdiff.min(),zdiff.min()))

		patch_areas = radiff.min()*decdiff.min()

		print('patching sky..')
		hist2d_c, hist2d_e = np.histogramdd(gal_coords[:,:2], bins=[redg,dedg])
		patch_idx = np.array(np.where(hist2d_c!=0)).T # identify populated patches

		patch_cuts = np.empty([patch_idx.shape[0], len(gal_coords)])
		pix_cuts = np.empty([patch_idx.shape[0], len(lpix_ra)])
		for c, (i,j) in enumerate(patch_idx):
			edges = redg[i], redg[i+1], dedg[j], dedg[j+1], 0., 0.5
			cut = betwixt(edges) # returns function cut(), which takes ra,dec,z and returns boolean array
			patch_cuts[c] = cut(ra,dec,z)
			pix_cuts[c] = cut(lpix_ra, lpix_dec, np.ones_like(lpix_ra)*0.1)

		npixLost = np.array([np.count_nonzero(i) for i in pix_cuts])
		pixar = hp.nside2pixarea(2048,degrees=True)
		pixperpatch = patch_areas/pixar
		patch_weights = 1-(npixLost/pixperpatch)
		print('patch weights (check none << 0) :\n',patch_weights)
		patch_weights = np.where(patch_weights>0,patch_weights,0)

		edges = (redg,dedg,zedg)
		return patch_cuts, patch_weights, patch_idx, edges

	def make_patches(self, fits_cat, patch_cuts):
		# apply patch cuts to fitsdata
		patches = np.array([fits_cat[np.array(i,dtype=bool)] for i in patch_cuts]) # array of populated patches, where each is a fits-table
		return patches

	def filter_patches(self, patches, patch_cuts, patch_weights, patch_idx, edges, lilbin_division=4, occupn_threshold=0.8):
		print('filtering survey-edges..')
		occ_fracs = np.empty(len(patches))
		pwei = np.empty(len(patches))
		occupied_area = 0
		retained_occ_area = 0
		area = lilbin_division**2 # = 16 (default) in arbitrary units
		norm_area = area*len(patches)

		for i, patch in enumerate(patches):
			# subdivide each patch for occupations
			coords = np.column_stack((patch[self.cols[0]], patch[self.cols[1]], patch[self.cols[2]]))
			r,d = patch_idx[i] # indices patch-edges
			patch_edges = (edges[0][r], edges[0][r+1], edges[1][d], edges[1][d+1])
			lilbins = (np.linspace(patch_edges[0], patch_edges[1], lilbin_division+1), np.linspace(patch_edges[2], patch_edges[3], lilbin_division+1))
                        fine_bins = (np.linspace(patch_edges[0], patch_edges[1], 11), np.linspace(patch_edges[2], patch_edges[3], 11))
			patch_hist2d = np.histogramdd(coords[:,:2], bins=lilbins)
			fine_hist2d = np.histogramdd(coords[:,:2], bins=fine_bins)

			occupn_frac = np.sum(patch_hist2d[0]!=0) / area
			pweight = np.sum(fine_hist2d[0]>50) / 100.
			occ_fracs[i] = occupn_frac
			pwei[i] = pweight**2
			occupied_area += pweight*area
			if occupn_frac>occupn_threshold:
				retained_occ_area += pweight*area

		if self.do_sdss:
			patch_filter = occ_fracs>=occupn_threshold
			
		else:
			patch_filter = np.ones_like(occ_fracs, dtype=bool)

		highfracs = (pwei**0.5)[patch_filter]
		discarded_fraction = sum(patch_filter)/len(patch_filter)
		print('occupation threshold (SDSS): ', occupn_threshold,
			'\ntotal (populated) patches: ', len(patch_filter),
			'\ndiscarded edge-patches (SDSS): ', sum(~patch_filter),
			'\nremaining: ', sum(patch_filter),
			'\nmin | max | mean | (step) in patch occupations: %.4f | %.4f | %.4f | (%.4f)'%(highfracs.min(), highfracs.max(), np.mean(highfracs), 1/area),
			)

		filtered_area_ratio = retained_occ_area/occupied_area
		scale_factor = filtered_area_ratio**-0.5
		if not self.do_sdss:
			print('GAMA: no discarded patches')
			scale_factor = 1.
		else:
			print('SCALE JACKKNIFE ERRORS DOWN BY ~%.3f TO ACCOUNT FOR LOST (EDGE) PATCHES...!!'%scale_factor)

		patch_cuts = patch_cuts[patch_filter]
		patch_weights = patch_weights[patch_filter]
		if self.do_sdss:
			patch_weights = pwei[patch_filter]
			print('SDSS patch weights from occupations:\n', patch_weights)
		return np.array(patch_cuts, dtype=bool), patch_weights, scale_factor

	def make_sample_patchcuts(self, patch_cuts, sample_cut):
		# apply sample cuts ONE AT A TIME to patch cuts
		new_patch_cuts = np.empty_like(patch_cuts)
		for i in range(len(patch_cuts)):
			new_patch_cuts[i] = patch_cuts[i]&sample_cut
		return new_patch_cuts

	def make_cubes(self, patches, patch_weights, zedges):
		# slice patches in redshift, creating jackknife cubes
		print('slicing patches into cubes..')
		cubes = []
		cube_weights = []
		if self.do_sdss:
			zedges = zedges[zedges <= 0.25]
		for i, patch in enumerate(patches):
			patch_z = patch[self.cols[2]]
			zcuts = [ (patch_z>=zedges[j]) & (patch_z<=zedges[j+1]) for j in range(len(zedges)-1) ]
			for zcut in zcuts:
				cubes.append(patch[zcut])
				cube_weights.append(patch_weights[i])

		cubes = np.array(cubes)
		cube_weights = np.array(cube_weights)
		print('c :', cubes.shape)
		print('cw :', cube_weights.shape)
		cubes = cubes.flatten()
		cube_weights = cube_weights.flatten()

		return cubes, cube_weights

	def find_patch_edges(self, patch_cuts):
		random_cutter = []
		for i, patch_cut in enumerate(patch_cuts):
			ra, dec, z = self.data[self.cols[0]], self.data[self.cols[1]], self.data[self.cols[2]]
			ra, dec, z = ra[patch_cut], dec[patch_cut], z[patch_cut]
			p_cut = betwixt((ra.min(), ra.max(), dec.min(), dec.max(), z.min(), z.max()))
			random_cutter.append(p_cut)

		random_cutter = np.array(random_cutter)

		return random_cutter

# sdss = 0

# cat = fits.open( ['./DEI_GAMAv20.fits', './RMandelbaum_catalogs/lss.14full0.rfgr.removegals.fits'] [sdss] )[1].data
# cols = [ ['RA_GAMA', 'DEC_GAMA', 'Z_TONRY'], ['ra', 'dec', 'z'] ] [sdss]
# ra,dec,z = (cat['%s'%i] for i in cols)
# coords = np.column_stack((ra,dec,z))

# ranges = [ [ (i, j, -3., 3., 0., 0.5) for i,j in [ (129.,141.), (174.,186.), (211.5,223.5) ] ], None ] [sdss]
# cubes, cube_idx, bins = define_cubes(coords, 3., 3., 0.06, ranges=ranges) # if making cubes smaller, then lower lilbin_division(integer) also - cube filter very sensitive to this!
# cubes = filter_cubes(cubes, cube_idx, bins, cols=cols, lilbin_division=4, occupn_threshold=0.9)

# cols2 = [ ['RA_GAMA', 'DEC_GAMA', 'Z_TONRY', 'absmag_g_1', 'absmag_i_1'], ['ra', 'dec', 'z', 'rf_g-r'] ] [sdss]
# out_cat = np.zeros([2, len(cols2)+1])

# for c, cube in enumerate(cubes):
# 	cube_cols = np.column_stack((cube['%s'%i] for i in cols2))
# 	cube_cols = np.column_stack((cube_cols, [c]*len(cube)))
# 	out_cat = np.append(out_cat, cube_cols, axis=0)

# cols2 = [ ['# RA_GAMA', 'DEC_GAMA', 'Z_TONRY', 'absmag_g_1', 'absmag_i_1'], ['# ra', 'dec', 'z', 'rf_g-r'] ] [sdss]
# out_cat = out_cat[2:]
# ascii.write(out_cat, './3DJK_TEST.asc', names=cols2+['cube id'], delimiter='\t')





