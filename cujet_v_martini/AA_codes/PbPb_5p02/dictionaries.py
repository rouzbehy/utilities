from numpy import sqrt
#multiplicity = {'00-05':1762,
#                '05-10':1380,
#                '10-15':1088,
#                '00-10':(1762+1380)/2.,
#                '15-20':855.3,
#                '00-20':(1762+1380+1088+855.3)/4.,
#                '10-20':(1088+855.3)/2,
#                '30-50':(392.9+294.5+216.4+155.5)/4}
multiplicity = {'00-05':1762-147,
                '05-10':1380-113,
                '10-15':1088-93.4,
                '00-10':(1762+1380)/2.,
                '15-20':855.3-80.8,
                '00-20':(1762+1380+1088+855.3)/4.,
                '10-20':(1088+855.3)/2,
                '30-50':(392.9+294.5+216.4+155.5)/4}

exag_factor = {'cujet':2000,'martini':2000}

markers = {'alice':'P','cms':'o','atlas1':'*','atlas2':'s'}

eloss_colours = {'cujet':"#7cae7a", 'martini':"#d64045", 'no jet-medium':'orange'}
#'atlas2':"#0092cc", 'cms':"#eca400"}
#, "#422040"}


channel_colors = {'jet-medium-cujet'   : "#7cae7a",
                  'jet-medium-martini' : "#d64045",
                  #'brem'   : "#d64045",
                  'prompt' : "#0092cc", 
                  'thermal': "#eca400", 
                  'preEq'  : "#422040"}

channel_linestyles = {'jet-medium-cujet': "solid",
                      'jet-medium-martini': "solid",
                      'prompt': "dashed", 
                      'thermal': "dotted", 
                      'preEq': "dashdot"}

JF_cent_labels = {"C0-5":"00-05","C20-30":"20-30","C30-40":"30-40"}

## Column names for the experimental results
## done to make plotting easier by using a uniform naming scheme
cols_CMS_chRAA = {'PT [GEV]':'x',
                      'PT [GEV] LOW':'xlow',
                      'PT [GEV] HIGH':'xhigh','RAA':'y','stat +':'dy_stat+',
                      'stat -':'dy_stat-',
                      'sys +':'dy_syst+',
                      'sys -':'dy_syst-',
                      'sys,TAA +':'dy_systTAA+',
                      'sys,TAA -':'dy_systTAA-','sys,lumi +':'dy_systLum+',
                      'sys,lumi -':'dy_systLum-'}
cols_ALICE_chRAA = {'PT [GEV]':'x','PT [GEV] LOW':'xlow','PT [GEV] HIGH':'xhigh',
                    'RAA':'y','stat +':'dy_stat+','stat -':'dy_stat+','sys +':'dy_syst+','sys -':'dy_syst-'}


## hadronic ratio observable:
## kaon to pions and protons to pions
colnames_kaon_to_pion={'$p_{T}$ [$GeV/c$]':'x',
                       '$p_{T}$ [$GeV/c$] LOW':'xlow',
                       '$p_{T}$ [$GeV/c$] HIGH':'xhigh',
                       '$(K^{+}+K^{-})/(\pi^{+}+\pi^{-})$':'y',
                       'stat. +':'dy_stat+',
                       'stat. -':'dy_stat-',
                       'syst. +':'dy_syst1+',
                       'syst. -':'dy_syst1-',
                       'syst. uncorr. +':'dy_syst2+',
                       'syst. uncorr. -':'dy_syst2-'}

colnames_proton_to_pion={'$p_{T}$ [$GeV/c$]':'x',
                       '$p_{T}$ [$GeV/c$] LOW':'xlow',
                       '$p_{T}$ [$GeV/c$] HIGH':'xhigh',
                       '$(p+pbar)/(\pi^{+}+\pi^{-})$':'y',
                       'stat. +':'dy_stat+',
                       'stat. -':'dy_stat-',
                       'syst. +':'dy_syst1+',
                       'syst. -':'dy_syst1-',
                       'syst. uncorr. +':'dy_syst2+',
                       'syst. uncorr. -':'dy_syst2-'}


pp_colnames_kaon_to_pion={'$p_{T}$ [$GeV/c$]':'x',
                       '$p_{T}$ [$GeV/c$] LOW':'xlow',
                       '$p_{T}$ [$GeV/c$] HIGH':'xhigh',
                       '$(K^{+}+K^{-})/(\pi^{+}+\pi^{-})$':'y',
                       'stat. +':'dy_stat+',
                       'stat. -':'dy_stat-',
                       'syst. +':'dy_syst+',
                       'syst. -':'dy_syst-'}

pp_colnames_proton_to_pion={'$p_{T}$ [$GeV/c$]':'x',
                       '$p_{T}$ [$GeV/c$] LOW':'xlow',
                       '$p_{T}$ [$GeV/c$] HIGH':'xhigh',
                       '$(p+pbar)/(\pi^{+}+\pi^{-})$':'y',
                       'stat. +':'dy_stat+',
                       'stat. -':'dy_stat-',
                       'syst. +':'dy_syst+',
                       'syst. -':'dy_syst-'}

