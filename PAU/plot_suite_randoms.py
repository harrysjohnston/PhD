from functions import *

rands = []
args = []
gama = paus = False
for rf in sys.argv[1:]:
	if 'K1000' in rf or 'gama' in rf:
		args.append(dict(cat='K1000_SML.fits', zspec='Z_TONRY', zphot='zphot_2', mag='PETROMAG_R', id='CATAID', zlim=0.6, red_seq='red_sequence'))
		gama = True
	elif 'PAU' in rf or 'pau' in rf:
		args.append(dict(cat='PAUS_KSB.fits', zspec='ZBEST', zphot='bcnz_zb', mag='mag_i', id='numeric_id', zlim=1.3, red_seq='red_sequence_Cigale_3cluster'))
		paus = True
	rands.append(rf)

if gama and paus:
	print 'best not to mix GAMA & PAUS'
	sys.exit()

#cat = fopen(args[0]['cat'])

#ipython.magic('run plot_zphot_vs_zspec_contours.py %s %s'%(' '.join([args[0]['cat'], args[0]['zspec'], args[0]['zphot'], args[0]['mag'], args[0]['id']]), ' '.join(rands)))
##ipython.magic('run plot_nofz.py %s %s'%(' '.join([args[0]['cat'], args[0]['zphot'], args[0]['mag'], args[0]['id']]), ' '.join(rands)))
##myhist(nn(cat[args[0]['zspec']]),lw=2,ls='--',color='k',range=(-0.01,args[0]['zlim']))
##for rf in rands:
##	ipython.magic('run plot_rb_nofz.py %s %s'%(' '.join([args[0]['cat'], args[0]['zphot'], args[0]['mag'], args[0]['id'], args[0]['red_seq']]), rf))
##	myhist(nn(cat[args[0]['zspec']][cat[args[0]['red_seq']]]),lw=2,ls='--',color='r',range=(-0.01,args[0]['zlim']))
##	myhist(nn(cat[args[0]['zspec']][~cat[args[0]['red_seq']]]),lw=2,ls='--',color='b',range=(-0.01,args[0]['zlim']))
#for i in range(len(rands)):
#	ipython.magic('run plot_nzspec_zphot.py %s %s %s'%(' '.join([args[i]['cat'],args[i]['zphot'],args[i]['zspec']]), rands[i], ' '.join([args[i]['id'],args[i]['mag']])))

#os.system('python plot_zphot_vs_zspec_contours.py %s %s'%(' '.join([args[0]['cat'], args[0]['zspec'], args[0]['zphot'], args[0]['mag'], args[0]['id']]), ' '.join(rands)))
for i in range(len(rands)):
	os.system('python plot_nzspec_zphot.py %s %s %s'%(' '.join([args[i]['cat'],args[i]['zphot'],args[i]['zspec']]), rands[i], ' '.join([args[i]['id'],args[i]['mag']])))

