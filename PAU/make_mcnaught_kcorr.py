# coding: utf-8
from functions import *
cat = fopen('../SMLambdarApMatchedPhotom.fits')
gr = cat['absmag_g'] - cat['absmag_r']
grbin = np.linspace(0,1,8)
grbin[0] = gr.min()
grbin[-1] = gr.max()
poly = pd.read_csv('mcnaught_poly', delimiter=' ')
zgrid = np.linspace(0, 0.6, 61)
kc = []
for i in range(7):
    kz = poly['a0'][i] * zgrid**4. + poly['a1'][i] * zgrid**3. + poly['a2'][i] * zgrid**2. + poly['a3'][i] * zgrid + poly['a4'][i]
    kz[0] = 0.
    kc.append(kz)
    median_gr = np.median(gr[(gr > grbin[i]) & (gr < grbin[i+1])])
    plt.plot(zgrid, kz, label=median_gr)
plt.legend()
plt.xlim(0,0.5)
plt.ylim(-0.2, 1.2)
plt.show()
plt.tight_layout()

colbins = np.digitize(gr, grbin, right=1) - 1
idcol, zcol, mr_col = [], [], []
for zi in range(len(zgrid)):
    idcol.append(cat['CATAID'])
    zcol.append([zgrid[zi]]*len(cat['CATAID']))
    kzi = np.array(kc)[colbins][:, zi]
    mri = 10**(-kzi/2.5)
    mr_col.append(mri)
idcol = np.array(idcol).flatten()
zcol = np.array(zcol).flatten()
mr_col = np.array(mr_col).flatten()
kcorr = np.column_stack((idcol, zcol, mr_col))  

with h5py.File('../SMLambdarApMatchedPhotom_McNaught.kcorrs', 'w') as f:
    f.create_dataset('ID', data=idcol)
    f.create_dataset('z', data=zcol)
    f.create_dataset('maggy_ratio', data=mr_col)
    f.close()
    
