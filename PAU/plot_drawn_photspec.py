# coding: utf-8
from functions import *
catpath = 'PAUS_KSB.fits'
cat = fopen(catpath)

f, ax = plt.subplots(1, 3, figsize=(16,5))

plt.sca(ax[0])
c1 = cat[~np.isnan(cat['ZBEST']) & (cat['ZBEST']!=-1)]
zhist2d = np.histogram2d(c1['ZBEST'], c1['bcnz_zb'], bins=50)
plt.imshow(log10(zhist2d[0].T), origin='lower', aspect='auto', cmap='bone_r',
    extent=[c1['ZBEST'].min(),c1['ZBEST'].max(),c1['bcnz_zb'].min(),c1['bcnz_zb'].max()])
plt.title('PAUSW3-DEEP2 matches')
plt.xlabel(r'$z_{\rm{spec.}}$')      
plt.ylabel(r'$z_{\rm{phot.}}$')

pdens = 1.*zhist2d[0] / zhist2d[0].sum()
SZ = np.repeat(midpoints(zhist2d[1]), 50).reshape(50,50)
PZ = np.repeat(midpoints(zhist2d[2]), 50).reshape(50,50).T
draw = np.random.choice(np.arange((len(zhist2d[1])-1)*(len(zhist2d[2])-1)), p=pdens.flatten(), size=len(cat), replace=True)
draw_SZ = SZ.flatten()[draw] + (np.random.rand(len(draw))-0.5)/15.
draw_PZ = PZ.flatten()[draw] + (np.random.rand(len(draw))-0.5)/15.
drawhist2d = np.histogram2d(draw_SZ, draw_PZ, bins=100)
plt.sca(ax[2])
plt.imshow(log10(drawhist2d[0].T), origin='lower', aspect='auto', cmap='bone_r',
    extent=[draw_SZ.min(),draw_SZ.max(),draw_PZ.min(),draw_PZ.max()])
plt.xlabel(r'$z_{\rm{spec.}}$')
plt.ylabel(r'$z_{\rm{phot.}}$')
plt.title('draws from P(spec, phot)')

plt.sca(ax[1])
plt.hist(c1['ZBEST'], alpha=0.4, bins='auto', density=1, label='spec-z')
plt.hist(c1['bcnz_zb'], alpha=0.4, bins='auto', density=1, label='bcnz2')
myhist(draw_SZ, color='C0', lw=1.6, label='spec-draw')
myhist(draw_PZ, color='C1', lw=1.6, label='phot-draw')
plt.legend()
plt.xlabel('$z$')
plt.tight_layout()
plt.show()

t = Table(cat)
t['pz_draw'] = draw_PZ
t['sz_draw'] = draw_SZ
t.write(catpath, overwrite=1)






