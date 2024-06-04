from matplotlib.colors import TABLEAU_COLORS as tab

colors = { 0: "#422040",
           1: "#d64045",
           2: "#0092cc", 
           3: "#eca400", 
           4: "#7cae7a"}

rate_colours = {'NP' : "#CC3311",
                'NLO': "#33BBEE",
                'LO' : "#EE7733"}
eloss_colours = {'cujet':"#7cae7a", 'martini':"#d64045"}
markers_experiment = {'alice':'P','cms':'o','atlas1':'*','atlas2':'s'}
markers_cone_size = {0:'s',1:'P',2:'d'}

pgrid_info = {"hspace": 0.043, 
              "left": 0.09, 
              "right": 0.979,
              "bottom": 0.09,
              "top": 0.95}

pgrid_info_2 = {"height_ratios": (3, 1), 
               "hspace": 0.043, 
               "left": 0.07, 
               "right": 0.98,
               "bottom": 0.09,
               "top": 0.95}

alpha_s_vals = {'LO':0.280,\
                'NLO':0.242,\
                'NP':0.260}

## Jets:
colnames_ATLAS_spec_pp = {'PT [GEV]' : 'x', 
                        'PT [GEV] LOW': 'xlow',
                        'PT [GEV] HIGH': 'xhigh',
                        'D2(SIG)/DPT/DYRAP [NB/GEV]': 'y',
                        'stat +': 'dy_stat+',
                        'stat -': 'dy_stat-',
                        'sys,uncorr +': 'dy_syst1+',
                        'sys,uncorr -': 'dy_syst1-',
                        'sys,corr +' : 'dy_syst+',
                        'sys,corr -' : 'dy_syst-',
                        'sys,luminosity uncertainty +': 'dy_syst3+',
                        'sys,luminosity uncertainty -': 'dy_syst3-'}

colnames_ALICE_spec_pp = {'PTjet [GEV]':'x',
                          'PTjet [GEV] LOW':'xlow',
                          'PTjet [GEV] HIGH':'xhigh',
                          '$d^{2}\sigma/dp_{T}d\eta$ [mb $c$/GeV]':'y',
                          'stat +':'dy_stat+',
                          'stat -':'dy_stat-',
                          'sys,correlated +':'dy_syst+',
                          'sys,correlated -':'dy_syst-',
                          'sys,shape +':'dy_shape+',
                          'sys,shape -':'dy_shape-',
                          'sys,lumi +':'dy_lumi+',
                          'sys,lumi -':'dy_lumi-'}

colnames_CMS_spec_pp = {'PT [GEV]':'x',
                        'PT [GEV] LOW':'xlow',
                        'PT [GEV] HIGH':'xhigh',
                        'd^2 sigma / d PT d ETA':'y',
                        'stat +':'dy_stat+',
                        'stat -':'dy_stat-',
                        'sys +':'dy_syst+',
                        'sys -':'dy_syst-',
                        'global +':'dy_global+',
                        'global -':'dy_global-'}


colnames_ALICE_spec = {'PT [GEV/C]':'x',
                       'PT [GEV/C] LOW':'xlow',
                       'PT [GEV/C] HIGH':'xhigh',
                       'D2N/DPT/DETA/NEVT/NCOLL [1/(GEV/C)]':'y',
                       'stat +': 'dy_stat+',
                       'stat -': 'dy_stat-',
                       'sys,shape +': 'dy_syst1+',
                       'sys,shape -': 'dy_syst1-',
                       'sys,corr +': 'dy_syst+',
                       'sys,corr -': 'dy_syst-',
                       'sys,$N_{\mathrm{coll}}$ uncertainty +':'dy_syst3+',
                       'sys,$N_{\mathrm{coll}}$ uncertainty -': 'dy_syst3-'}
colnames_CMS_spec = {'Jet $p_T$ [GEV/c]':'x',
                     'Jet $p_T$ [GEV/c] LOW':'xlow',
                     'Jet $p_T$ [GEV/c] HIGH':'xhigh',
                     '$10^2 * 1 /(TAA*Nevt) * d^2 N / dp_T d\eta$ [nb/GEV/c]':'y',
                     'stat +':'dy_stat+',
                     'stat -':'dy_stat-',
                     'sys +':'dy_syst+',
                     'sys -':'dy_syst-'}

colnames_CMS_RAA = {'PT [GEV]':'x',
                    'PT [GEV] LOW':'xlow',
                    'PT [GEV] HIGH':'xhigh',
                    'RAA':'y',
                    'stat +':'dy_stat+',
                    'stat -': 'dy_stat-',
                    'sys +':'dy_syst+',
                    'sys -': 'dy_syst-',
                    'TAA +':'dy_syst1+',
                    'TAA -':'dy_syst1-',
                    'lumi +':'dy_syst2+',
                    'lumi -':'dy_syst2-'}

colnames_ATLAS_RAA = {"PT [GEV] [GeV]":'x',
                      "PT [GEV] [GeV] LOW":'xlow',
                      "PT [GEV] [GeV] HIGH": 'xhigh',
                      'RAA':'y',
                      "stat +":'dy_stat+',
                      "stat -":'dy_stat-',
                      "sys +":'dy_syst+',
                      "sys -":'dy_syst-',
                      "sys,luminosity uncertainty +":'dy_syst2+',
                      "sys,luminosity uncertainty -":'dy_syst2-',
                      "sys,TAA uncertainty +"       :'dy_syst3+',
                      "sys,TAA uncertainty -"       :'dy_syst3-'}

colnames_ATLAS_FF_pT = {'PT [GEV]':'x',
                	    'PT [GEV] LOW':'xlow',
                        'PT [GEV] HIGH':'xhigh',
                        'DN/DPT (P P)':'y',
                        'stat +':'dy_stat+',
                        'stat -':'dy_stat-',
                        'sys +':'dy_syst+',
                        'sys -':'dy_syst-'}
colnames_ATLAS_FF_z = { 'Z':'x',
                        'Z LOW':'xlow',
                        'Z HIGH':'xhigh',
                        'DN/DZ (P P)':'y',
                        'stat +':'dy_stat+',
                        'stat -':'dy_stat-',
                        'sys +': 'dy_syst+',
                        'sys -': 'dy_syst-' }
