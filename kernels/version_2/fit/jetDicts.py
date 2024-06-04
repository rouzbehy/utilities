from matplotlib.colors import TABLEAU_COLORS as tab

colors = { 0 : tab['tab:blue'],
           1 : tab['tab:orange'],
           2 : tab['tab:green'],
           3 : tab['tab:red'],
           4 : tab['tab:purple']}

rate_colours = {'NP' : "#CC3311",
                'NLO': "#33BBEE",
                'LO' : "#EE7733"}


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
                        'sys,corr +' : 'dy_syst2+',
                        'sys,corr -' : 'dy_syst2-',
                        'sys,luminosity uncertainty +': 'dy_syst3+',
                        'sys,luminosity uncertainty -': 'dy_syst3-'}
colnames_ALICE_spec_pp = {'PT [GEV/C]' : 'x',
                        'PT [GEV/C] LOW' : 'xlow',
                        'PT [GEV/C] HIGH': 'xhigh',
                        'D2N/DPT/DETA/NEVT [1/(GEV/C)]' : 'y',
                        'stat +' : 'dy_stat+',
                        'stat -' : 'dy_stat-',
                        'sys,uncorr +': 'dy_syst1+',
                        'sys,uncorr -': 'dy_syst1-',
                        'sys,total inelastic cross section uncertainty +' : 'dy_syst2+',
                        'sys,total inelastic cross section uncertainty -' : 'dy_syst2-'}
colnames_CMS_spec_pp = {'Jet $p_T$ [GEV/c]': 'x',
                      'Jet $p_T$ [GEV/c] LOW': 'xlow',
                      'Jet $p_T$ [GEV/c] HIGH': 'xhigh',
                      '$d^2 \sigma / dp_T d\eta$ [nb/GEV/c]': 'y',
                      'stat +': 'dy_stat+',
                      'stat -': 'dy_stat-',
                      'sys +': 'dy_syst+',
                      'sys -': 'dy_syst-'}

colnames_ALICE_spec = {'PT [GEV/C]':'x',
                       'PT [GEV/C] LOW':'xlow',
                       'PT [GEV/C] HIGH':'xhigh',
                       'D2N/DPT/DETA/NEVT/NCOLL [1/(GEV/C)]':'y',
                       'stat +': 'dy_stat+',
                       'stat -': 'dy_stat-',
                       'sys,shape +': 'dy_syst1+',
                       'sys,shape -': 'dy_syst1-',
                       'sys,corr +': 'dy_syst2+',
                       'sys,corr -': 'dy_syst2-',
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


colnames_CMS_RAA = {'Jet $p_T$ [GEV/c]':'x',
                    'Jet $p_T$ [GEV/c] LOW':'xlow',
                    'Jet $p_T$ [GEV/c] HIGH':'xhigh',
                    'RAA':'y',
                    'stat +':'dy_stat+',
                    'stat -': 'dy_stat-',
                    'sys +':'dy_syst+',
                    'sys -': 'dy_syst-'}

colnames_ATLAS_RAA = {"PT [GEV]":'x',
                      "PT [GEV] LOW":'xlow',
                      "PT [GEV] HIGH": 'xhigh',
                      'RAA':'y',
                      "stat +":'dy_stat+',
                      "stat -":'dy_stat-',
                      "sys,corr +"                 :'dy_syst1+',
                      "sys,corr -"                  :'dy_syst1-',
                      "sys,uncorr +"                :'dy_syst2+',
                      "sys,uncorr -"                :'dy_syst2-',
                      "sys,luminosity uncertainty +":'dy_syst3+',
                      "sys,luminosity uncertainty -":'dy_syst3-',
                      "sys,TAA uncertainty +"       :'dy_syst4+',
                      "sys,TAA uncertainty -"       :'dy_syst4-'}
