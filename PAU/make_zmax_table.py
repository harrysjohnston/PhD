# coding: utf-8
from functions import *
from clone_randoms import *
rup = lambda x: int(np.ceil(x))

catpath = sys.argv[1]
zph_col = sys.argv[2]
zsp_col = sys.argv[3]
id_col = sys.argv[4]
Ndraw = int(sys.argv[5])
kcorr_file = sys.argv[6]
mlim = float(sys.argv[7])
mag_col = sys.argv[8]
Q = float(sys.argv[9])
zlimit = float(sys.argv[10])
dz = float(sys.argv[11])
mode = sys.argv[12] # bins, pergal or plot

cat = fopen(catpath)
cat = cat[np.argsort(cat[id_col])]
zph = cat[zph_col]
zsp = cat[zsp_col]
mag = cat[mag_col]

if mode == 'pergal' or mode == 'bins':
	plot_only = 0
if mode == 'plot':
	plot_only = 1
	mode = 'bins'

bins = np.linspace(0, zlimit, rup(zlimit/dz)+1)
nbin = len(bins)-1

if plot_only:
	f, ax = plt.subplots(rup(nbin**0.5), rup(nbin**0.5), sharex=True, sharey=True, figsize=(16,9))
	ai = iter(ax.flatten())

#for n in range(Ndraw):
def draw_zmax(n, mag=mag, bins=bins):
	zsp_draw_lookup = []

	for i in range(len(bins)-1):
		zph_cut = (zph > bins[i]) & (zph <= bins[i+1])
		zph_i = zph[zph_cut]
		zsp_i = zsp[zph_cut]
		catid_i = cat[id_col][zph_cut]

		hist_i, bin_i = np.histogram(nn(zsp_i), bins=bins)
		if all(hist_i == 0):
			continue

		norm_hist = 1.*hist_i / hist_i.sum()
		d_grid = midpoints(get_d(bin_i))

		if plot_only:
			plt.sca(next(ai))
			plt.plot(get_z(d_grid), norm_hist, label='zspec. dist')
			plt.axvline(bins[i], label='zphot. bin', c='k')
			plt.axvline(bins[i+1], c='k')
		
		zsp_draw = np.random.choice(d_grid, p=norm_hist, size=zph_cut.sum()) + \
					(np.random.rand(zph_cut.sum()) - 0.5) * np.diff(d_grid)[0]
		zsp_draw_lookup.append(np.column_stack((catid_i, zsp_draw)))

		if plot_only:
			hist_zspdraw, bin_zspdraw = np.histogram(get_z(nn(zsp_draw)), bins=bins)
			plt.plot(midpoints(bin_zspdraw), 1.*hist_zspdraw/hist_zspdraw.sum(), label='zspec. draw')
			plt.xlim(-0.02, bins.max()+0.02)

	if plot_only:
		plt.tight_layout()
		plt.subplots_adjust(wspace=0.03, hspace=0.05)
		ax.flatten()[0].legend(loc='best')
		plt.show()
		sys.exit()

	zsp_draw_lookup = np.concatenate(np.array(zsp_draw_lookup))
	zsp_draw_lookup = zsp_draw_lookup[np.argsort(zsp_draw_lookup[:, 0])]
	catid = zsp_draw_lookup[:, 0]
	zdraw = get_z(zsp_draw_lookup[:, 1])

	k_z, z_grid, mask = get_k_z(kcorr_file, catid)
	idcut = np.isin(cat[id_col], catid)
	if len(mag) > idcut.sum():
		mag = mag[idcut]

	max_redshift = fit_zmax(mlim, mag, zdraw, z_grid, k_z, Q=Q)
	max_redshift[mask] = -99.
	return np.column_stack((catid, max_redshift, zdraw))

def draw_zsp_pergal(i, mag=mag, bins=bins):

	zph_cut = (zph > zph[i]-dz) & (zph <= zph[i]+dz)
	zph_i = zph[zph_cut]
	zsp_i = zsp[zph_cut]
	catid_i = np.repeat(cat[id_col][i], Ndraw)

	hist_i, bin_i = np.histogram(nn(zsp_i), bins=bins)
	if all(hist_i == 0):
		return np.column_stack((catid_i, [-99.]*Ndraw))

	norm_hist = 1.*hist_i / hist_i.sum()
	d_grid = midpoints(get_d(bin_i))
	
	zsp_draw = np.random.choice(d_grid, p=norm_hist, size=Ndraw) + \
				(np.random.rand(Ndraw) - 0.5) * np.diff(d_grid)[0]
	zsp_draw = get_z(zsp_draw)
	zsp_draw_lookup = np.column_stack((catid_i, zsp_draw))
	return zsp_draw_lookup

def fit_zmax_pergal(n, zdraw, catid, mag=mag):
	sort = np.argsort(catid)
	k_z, z_grid, mask = get_k_z(kcorr_file, catid[sort])
	idcut = np.isin(cat[id_col], catid[sort])
	if len(mag) > idcut.sum():
		mag = mag[idcut]
	max_redshift = fit_zmax(mlim, mag, zdraw[sort], z_grid, k_z, Q=Q)
	max_redshift[mask] = -99.
	return np.column_stack((catid[sort], max_redshift, zdraw[sort]))

if plot_only:
	draw_zmax(0)
else:
	import multiprocessing as mp
	pool = mp.Pool(mp.cpu_count())
	if mode == 'bins':
		max_redshift_draws = np.array(pool.map(draw_zmax, range(Ndraw)))
	elif mode == 'pergal':
		zsp_draws = np.array(pool.map(draw_zsp_pergal, range(len(zph))))
		max_redshift_draws = np.array([pool.apply(fit_zmax_pergal,
										args=(n, zsp_draws[:,n::Ndraw,1].squeeze(), zsp_draws[:,n::Ndraw,0].squeeze()))
										for n in range(Ndraw)])
		true_max_redshift = np.array(pool.apply(fit_zmax_pergal, args=(0, np.where(zsp > 0, zsp, np.nan), zsp_draws[:,0::Ndraw,0].squeeze())))
	pool.close()

	catids = max_redshift_draws[0,:,0]
	zmax_table = max_redshift_draws[:,:,1]
	zspec_table = max_redshift_draws[:,:,2]

	if mode == 'bins':
		outfile = catpath.replace('.fits', '.dz%s.zmaxtable'%dz)
	elif mode == 'pergal':
		outfile = catpath.replace('.fits', '.pergaldz%s.zmaxtable'%dz)
	with h5py.File(outfile, 'w') as f:
		f.create_dataset('ID', data=catids)
		f.create_dataset('zmax', data=zmax_table)
		f.create_dataset('zspec', data=zspec_table)
		f.create_dataset('true_zmax', data=true_max_redshift)
		f.close()



