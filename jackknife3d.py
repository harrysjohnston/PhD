# coding: utf-8
from __future__ import print_function, division ; from astropy.io import fits, ascii ; import numpy as np ; import matplotlib.pylab as plt
from astropy.cosmology import Planck13
import healpy as hp
import gc

def betwixt((re1,re2,de1,de2,ze1,ze2)):
	def make_cut(ra,dec,z):
		return (ra>=re1)&(ra<=re2)&(dec>=de1)&(dec<=de2)&(z>=ze1)&(z<=ze2)
	return make_cut

def resample_data(fitsdata, sample_cuts, patchside=6, zcut=0.26, do_sdss=0, do_3d=1, cube_zdepth=0.06, largePi=0, bitmaskCut=None, occ_thresh=0.5, mask_path='/share/splinter/hj/PhD/pixel_weights.fits'):
	resampler = resampleTools(fitsdata, patchside, zcut, do_sdss, do_3d, sample_cuts, cube_zdepth=cube_zdepth, largePi=largePi)

	# identify masked pixel coordinates
	# lostpix_coords = resampler.find_lostpixels(bitmaskCut)
	# lostpix_coords = (np.array([1e3]),np.array([1e3]))

	if not do_sdss:
		# read pixel coordinates & weights
		pixel_coords = resampler.read_pixel_weights(mask_path=mask_path)
	else:
		pixel_coords = None

	# define patches & their degree of masking
	patch_cuts, patch_weights, patch_idx, edges = resampler.define_edgecuts(pixel_coords=pixel_coords)

	# apply patch cuts to data
	patches = resampler.make_patches(fitsdata, patch_cuts)

	# filter patches straddling survey-edges
	patch_cuts, patch_weights, error_scaling = resampler.filter_patches(patches, patch_cuts, patch_weights, patch_idx, edges,occupn_threshold=occ_thresh)

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
		# slice patches in z, creating jackknife cubes & random cube-cuts
		cubeData = []
		for i, spatches in enumerate(patchData):
			if i==0:
				cubes, cube_weights, random_cutter = resampler.make_cubes(spatches, patch_weights, edges[2], random_cutter=random_cutter)
			else:
				cubes, cube_weights = resampler.make_cubes(spatches, patch_weights, edges[2])
			cubeData.append(cubes)
	else:
		patchData = np.array(patchData)
		cubeData = patchData.copy()
		cube_weights = patch_weights.copy()
		del patchData, patch_weights

	cubeData = np.array(cubeData)
#	empty_cut = np.ones(len(cubeData), dtype=bool)
#	# rough cut against partially empty patches/cubes
#	for i, cube

	print('\nN cubes: ', len(cubeData[0]),
			'\nmin | max cube weights: %.4f | %.4f'%(cube_weights.min(), cube_weights.max()))

	return cubeData, cube_weights, error_scaling, random_cutter

class resampleTools:
	def __init__(self, fitsdata, patchside, zcut, do_sdss, do_3d, sample_cuts, cube_zdepth=0.06, largePi=0):
		self.cols = [ ['RA_GAMA', 'DEC_GAMA', 'Z_TONRY'], ['ra', 'dec', 'z'] ] [do_sdss]
		self.ranges = [ [ (i, j, -3., 3., 0., 0.6) for i,j in [ (129.,141.), (174.,186.), (211.5,223.5) ] ], None ] [do_sdss]
		self.data = fitsdata
		self.zcut = zcut
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
			#lostpixIDs = np.where(kidsBitmap==0, 0, 1)
			#lostpixIDs = np.arange(len(kidsBitmap))[lostpixIDs]
			bitmask_cut = [True]*len(kidsBitmap)

			print("PATCHES: CUTTING MASK!=0")
			bitmask_cut = np.where(kidsBitmap==0,True,False)

			if type(bitmask_[0]) != type(None):
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

	def read_pixel_weights(self, mask_path):
		# make 2d-array of pixel coords and weights, for patch_weighting
		mask_map = hp.read_map(mask_path)
		pix = np.arange(len(mask_map))
		ra, dec = hp.pix2ang(2048, pix, lonlat=True)
		z = np.ones_like(ra) * 0.2

		pixel_coords = np.column_stack(( ra, dec, z, mask_map ))
		return pixel_coords

	def define_edgecuts(self, pixel_coords=None): # if gama, pixel_coords has 4 cols = (ra, dec, z(dummy), weight) and npix-rows
		# given desired patch/box-sizes, divide sky into patches
		# return sets of patch-cuts for application to catalogs, with weights due to lost pixels

		gal_coords = np.column_stack((self.data[self.cols[0]], self.data[self.cols[1]], self.data[self.cols[2]]))
		ra, dec, z = gal_coords.T
		gama = not self.do_sdss

		print('defining cubes..')
		ra_num, dec_num, z_num = 360//self.ra_side +1, 180//self.dec_side +1, 0.6//self.z_side +1
		if self.largePi:
			z_num = 0.6//0.1 +1
		redg,dedg,zedg = np.linspace(0,360,ra_num), np.linspace(-90,90,dec_num), np.linspace(0,0.6,z_num)
		if type(self.ranges)!=type(None):
			print('GAMA equatorial: jackknifing specified ranges...')
			redg = []
			dec_num = 6.//self.dec_side +1
			dedg = np.linspace(-3., 3., dec_num)
			zedg = np.linspace(0., 0.6, z_num)
			if self.zcut != None:
				z_push = self.zcut - min(zedg, key=lambda x:abs(x - self.zcut))
				zedg = [zz for zz in (zedg + z_push) if (zz > 0.) & (zz < 0.6)]
			for edges in self.ranges:
				ra_num = (edges[1] - edges[0])//self.ra_side + 1
				redg += list(np.linspace(edges[0], edges[1], ra_num))

		redg,dedg,zedg = map(lambda x: np.array(x), [redg,dedg,zedg])
		if self.do_sdss:
			print('SDSS: forcing cube z-edges to 0., 0.125, 0.25')
			zedg = np.linspace(0., 0.5, 5)
		radiff, decdiff, zdiff = (np.diff(i) for i in [redg,dedg,zedg])
		print('ra | dec | z sides: %.2f | %.2f | %.3f'%(radiff.min(),decdiff.min(),zdiff.min()))
		self.ra_side, self.dec_side = radiff.min(), decdiff.min()

		# identify populated patches
		print('patching sky..')
		hist2d_c, hist2d_e = np.histogramdd(gal_coords[:,:2], bins=[redg,dedg])
		patch_idx = np.array(np.where(hist2d_c!=0)).T

		patch_cuts = np.empty( [patch_idx.shape[0], len(gal_coords)] )
		if gama:
			patch_pixel_weights = np.zeros( [patch_idx.shape[0], len(pixel_coords)] )
			patch_pixel_count = np.empty( patch_idx.shape[0] )
		for c, (i,j) in enumerate(patch_idx):
			edges = redg[i], redg[i+1], dedg[j], dedg[j+1], 0., 0.6
			cut = betwixt(edges) # returns function cut(), which takes ra,dec,z and returns boolean array
			patch_cuts[c] = cut(ra,dec,z)
			if gama:
				patch_pixel_weights[c] += np.where( cut(*pixel_coords.T[:3]), pixel_coords.T[3], 0 )
				patch_pixel_count[c] = np.sum( cut(*pixel_coords.T[:3]) )

		if gama:
			patch_weights = np.zeros( len(patch_pixel_weights) )
			for i in range(len(patch_weights)):
				patch_weights[i] = np.sum(patch_pixel_weights[i]) / np.float32( patch_pixel_count[i] )
		else:
			patch_weights = np.ones( patch_idx.shape[0] )
		print('===============\t CHECK THIS \t================\npatch weights:\n',patch_weights)

		edges = (redg,dedg,zedg)
		return patch_cuts, patch_weights, patch_idx, edges

	def make_patches(self, fits_cat, patch_cuts):
		# apply patch cuts to fitsdata
		patches = np.array([fits_cat[np.array(i,dtype=bool)] for i in patch_cuts]) # array of populated patches, where each is a fits-table
		return patches

	def filter_patches(self, patches, patch_cuts, patch_weights, patch_idx, edges, occupn_threshold=0.5):
		print('filtering survey-edges..')
		pwei = np.empty(len(patches))
		parea = np.empty(len(patches))
		survey_area = 0
		jackknife_area = 0
		Area = lambda a1,a2,d1,d2: (a2-a1)*(np.sin(d2)-np.sin(d1))

		for i, patch in enumerate(patches):
			# subdivide each patch for occupations
			coords = np.column_stack((patch[self.cols[0]], patch[self.cols[1]], patch[self.cols[2]]))
			r,d = patch_idx[i] # indices patch-edges
			patch_edges = (edges[0][r], edges[0][r+1], edges[1][d], edges[1][d+1])
			patch_edges_rad = [np.deg2rad(pe) for pe in patch_edges]
			fine_bins = (np.linspace(patch_edges[0], patch_edges[1], 11), np.linspace(patch_edges[2], patch_edges[3], 11))
			fine_hist2d = np.histogramdd(coords[:,:2], bins=fine_bins)

			pwei[i] = np.sum(fine_hist2d[0]!=0) / 100.
			parea[i] = Area(*patch_edges_rad) * (180./np.pi)**2
			survey_area += pwei[i] * parea[i]
			if pwei[i] >= occupn_threshold:
				jackknife_area += survey_area

		if self.do_sdss:
			patch_filter = pwei >= occupn_threshold
		else:
			patch_filter = np.ones_like(pwei, dtype=bool)

		highfracs = pwei[patch_filter]
		print('occupation threshold (SDSS): ', occupn_threshold,
			'\ntotal (populated) patches: ', len(patch_filter),
			'\ndiscarded edge-patches (SDSS): ', sum(~patch_filter),
			'\nremaining: ', sum(patch_filter),
			'\nmin | max | mean | (step) in patch occupations: %.4f | %.4f | %.4f | (%.4f)'%(highfracs.min(), highfracs.max(), np.mean(highfracs), 1/100.),
			)

		area_scaling = jackknife_area/survey_area
		scale_factor = area_scaling**-0.5
		if not self.do_sdss:
			print('GAMA: no discarded patches')
			scale_factor = 1.
		else:
			print('retained vs. occupied area fraction: %.3f'%area_scaling)
			print('SCALE JACKKNIFE ERRORS DOWN BY ~%.3f TO ACCOUNT FOR LOST (EDGE) PATCHES...!!'%scale_factor)

		patch_cuts = patch_cuts[patch_filter]
		patch_weights = patch_weights[patch_filter]
		if self.do_sdss:
			patch_weights = pwei[patch_filter] * parea[patch_filter] / (self.ra_side * self.dec_side)
			print('SDSS patch weights from area * occupation:\n', patch_weights)
		return np.array(patch_cuts, dtype=bool), patch_weights, scale_factor

	def make_sample_patchcuts(self, patch_cuts, sample_cut):
		# apply sample cuts ONE AT A TIME to patch cuts
		new_patch_cuts = np.empty_like(patch_cuts)
		for i in range(len(patch_cuts)):
			new_patch_cuts[i] = patch_cuts[i] & sample_cut
		return new_patch_cuts

	def make_cubes(self, patches, patch_weights, zedges, random_cutter=None):
		# slice patches in redshift, creating jackknife cubes
		print('slicing patches into cubes..')
		cubes = []
		cube_weights = []
		new_random_cutter = []
		if self.do_sdss:
			zedges = zedges[zedges <= 0.25]
		for i, patch in enumerate(patches):
			patch_z = patch[self.cols[2]]
			zcuts = [ (patch_z>=zedges[j]) & (patch_z<=zedges[j+1]) for j in range(len(zedges)-1) ]
			for k, zcut in enumerate(zcuts):
                                cubes.append( patch[zcut] )
                                cube_weights.append( patch_weights[i] )
				if type(random_cutter)!=type(None):
					# make cut() function for random cubing - need to discard cubes
					cube_cut = betwixt((-np.inf, np.inf, -np.inf, np.inf, zedges[k], zedges[k+1]))
					new_random_cutter.append( (random_cutter[i], cube_cut) )

		cubes = np.array(cubes)
		cube_weights = np.array(cube_weights)
		cubes = cubes.flatten()
		cube_weights = cube_weights.flatten()

		new_random_cutter = np.array(new_random_cutter)

		if type(random_cutter)!=type(None):
			return cubes, cube_weights, new_random_cutter
		else:
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


