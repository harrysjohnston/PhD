from functions import *
from tqdm import tqdm

cat = fopen('PAUS_KSB.fits')
ra = cat['alpha_j2000']
dec = cat['delta_j2000']

rra = []
rdec = []
for i in tqdm(range(50),ncols=100):
	rra.append(ra + (np.random.rand(len(ra))-0.5)*(1./60.))
	rdec.append(dec + (np.random.rand(len(dec))-0.5)*(1./60.))

rra = np.array(rra).flatten()
rdec = np.array(rdec).flatten()

c = np.random.rand(len(rra)) < 1
c1 = np.random.rand(len(ra)) < 1
plt.figure(figsize=(17,9))
plt.scatter(rra[c], rdec[c], s=0.1)
plt.scatter(ra[c1], dec[c1], s=0.3)
plt.show()

mask = radec_to_map(rra,rdec,nside=2048)!=0
hp.mollzoom(mask)
hp.write_map('W3_2048mask2.fits', mask, overwrite=1)


