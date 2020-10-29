dd = read_files('wgp*dat', dire='OUTPUTS_PAUS_KSBarea_fg1fg2_PW1/', asc=1)                                                                               
     ...: rp = dd['wgp_red_ra10.dat']['rnom']                                                                       
     ...: area_perc = [20,30,50,80]                                                                                                                                
     ...: cm = plt.cm.viridis(np.linspace(0,1,len(area_perc)))                                                                                                     
     ...: del e1                                                                                             
     ...: eratio = []                                                                                             
     ...: for i, k in enumerate(['wgp_red_ra%s.dat'%a for a in area_perc]):                                         
     ...:     r, e = dd[k]['rnom'], dd[k]['wgplus_jackknife_err']
     ...:     plt.plot(r, e, c=cm[i], label='%s%% of area'%area_perc[i])
     ...:     try:                               
     ...:         eratio.append(e/e1)                                                                 
     ...:         print list(np.round(e/e1, 3))              
     ...:     except:            
     ...:         pass                          
     ...:     e1 = e.copy()                           
     ...: plt.xscale('log')                                  
     ...: plt.yscale('log')                                                   
     ...: plt.ylabel(r'$\sigma(w_{\rm{g+}})$')  
     ...: plt.xlabel(r'$r_{p}$')
     ...: plt.legend()                                       
     ...: plt.tight_layout()  
