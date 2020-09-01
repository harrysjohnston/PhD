from functions import *

if len(sys.argv) < 9:
	print '\n\
Must give:\n\
1: zmaxtable.h5 path\n\
2: catalog path\n\
3-7: column names for photo-z, detection-mag, ID, RA and DEC (5 args)\n\
8: index of zspec/zmax table to draw\n\
9: output file identifier\n'
	sys.exit()

zmaxtable_h5 = h5py.File(sys.argv[1], 'r')
zmaxtable = zmaxtable_h5['zmax'][:]
zspectable = zmaxtable_h5['zspec'][:]
zmax_id = zmaxtable_h5['ID'][:]
zspectable = zspectable[:, np.argsort(zmax_id)]
zmaxtable = zmaxtable[:, np.argsort(zmax_id)]
zmaxtable_h5.close()

cat = fopen(sys.argv[2])
#cols = ['bcnz_zb', 'mag_i', 'numeric_id', 'alpha_j2000', 'delta_j2000']
# photo-z detection-mag unique-ID RA DEC
cols = sys.argv[3:8]
t = Table(dict(zip(cols, [cat[col] for col in cols])))
t = t[np.argsort(t[cols[2]])]

i = int(sys.argv[8])
outid = sys.argv[9]
t['zspec'] = zspectable[i]
t['zmax'] = zmaxtable[i]
t.write(outid+'%s.fits'%str(i).zfill(3),overwrite=1)

