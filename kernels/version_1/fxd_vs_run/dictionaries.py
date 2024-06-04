colnames_CMS_spectra_had = {'PT [GEV]' : 'x',
                            'PT [GEV] LOW': 'xlow',
                            'PT [GEV] HIGH': 'xhigh',
                            'E*D3(N)/DP**3 [GEV**-2]' :'y',
                            'stat +' : 'dy_stat+',
                            'stat -' : 'dy_stat-',
                            'sys +' : 'dy_syst+',
                            'sys -': 'dy_syst-'}
colnames_ALICE_spectra_had = {'PT [GEV]' : 'x', 
                              'PT [GEV] LOW': 'xlow',
                              'PT [GEV] HIGH': 'xhigh',
                              '(1/Nevt)*(1/*(2PI*PT))*D2(N)/DETARAP/DPT [GEV**-2]': 'y',
                              'stat +': 'dy_stat+',
                              'stat -': 'dy_stat-',
                              'sys +': 'dy_syst+',
                              'sys -': 'dy_syst-'}
colnames_ATLAS_spectra_had = {'PT [GeV]' : 'x',
                              'PT [GeV] LOW' : 'xlow',
                              'PT [GeV] HIGH': 'xhigh',
                              '1/NEVT 1/(2*PI*PT) 1/<TAA> D2(N)/DPT/DETA [MB/GeV2]' : 'y',
                              'stat +' : 'dy_stat+',
                              'stat -' : 'dy_stat-',
                              'sys +' : 'dy_syst+',
                              'sys -' : 'dy_syst-'}

name_dicts = {'atlas1':colnames_ATLAS_spectra_had,
              'atlas2':colnames_ATLAS_spectra_had,
              'alice' :colnames_ALICE_spectra_had,
              'cms'   :colnames_CMS_spectra_had}

spectra_fname_ATLAS1 = {'00-05':'../../experiment/PbPb_2p76/charged/ATLAS/ATLAS_ch_had_spec_eta_1_PbPb_00_05.csv'}

spectra_fname_ATLAS2 = {'00-05':'../../experiment/PbPb_2p76/charged/ATLAS/Table2.csv',
                        '05-10':'../../experiment/PbPb_2p76/charged/ATLAS/Table3.csv',
                        '10-20':'../../experiment/PbPb_2p76/charged/ATLAS/Table4.csv'}

spectra_fname_ALICE = {'00-05': '../../experiment/PbPb_2p76/charged/ALICE/Table1.csv',
                       '05-10': '../../experiment/PbPb_2p76/charged/ALICE/Table2.csv',
                       '10-20': '../../experiment/PbPb_2p76/charged/ALICE/Table3.csv',
                       '00-20': '../../experiment/PbPb_2p76/charged/ALICE/Table4.csv'}

spectra_fname_CMS = {'00-05': "../../experiment/PbPb_2p76/charged/CMS/Table3_1.csv",
                     '05-10': "../../experiment/PbPb_2p76/charged/CMS/Table3_2.csv"}

spectrum_dicts = {'alice': spectra_fname_ALICE,
                  'atlas1':spectra_fname_ATLAS1,
                  'atlas2':spectra_fname_ATLAS2,
                  'cms':spectra_fname_CMS}

etadict = {'alice':0.8,
           'atlas2':2.0,
           'atlas1':1.0,
           'cms':1.0}

## hydro files:
hydro_spec_fname = {'00-05' : "../hydro_averages/PbPb_2p76/charged_hadron_differential_observables_ALICE_00-05.dat",
                    '05-10' : "../hydro_averages/PbPb_2p76/charged_hadron_differential_observables_ALICE_05-10.dat",
                    '10-20' : "../hydro_averages/PbPb_2p76/charged_hadron_differential_observables_ALICE_10-20.dat"}
hydro_column_names = {'pT' : 'pT',
                      'dN/(2pi dy pT dpT)' : 'Nch',
                      'dN/(2pi dy pT dpT)_err' : 'dNch'}

## proton-proton baseline:
pp_fname = {'alice' : "../pp_baseline_calc/PbPb_2p76/new_pp_baseline/alice_compiled_charged_hadron.csv",
            'atlas2': "../pp_baseline_calc/PbPb_2p76/new_pp_baseline/atlas2_compiled_charged_hadron.csv",
            'atlas1': "../pp_baseline_calc/PbPb_2p76/new_pp_baseline/atlas1_compiled_charged_hadron.csv",
            'cms'   : "../pp_baseline_calc/PbPb_2p76/new_pp_baseline/cms_compiled_charged_hadron.csv"}
colnames_pp = {'pT'         : 'pT',
               'dpT'        : 'dpT',
               'charged'    : 'Nch',
               'chargederr' : 'dNch'}
## RAA results:
raa_fname_atlas1 = {'00-05':"../../experiment/PbPb_2p76/charged/ATLAS/Table49.csv",## change this to 49 for eta < 1 instead of 2
}

raa_fname_atlas2 = {'00-05':"../../experiment/PbPb_2p76/charged/ATLAS/Table33.csv",## change this to 49 for eta < 1 instead of 2
                    '05-10':"../../experiment/PbPb_2p76/charged/ATLAS/Table34.csv",
                    '10-20':"../../experiment/PbPb_2p76/charged/ATLAS/Table35.csv"}

raa_fname_alice = {'00-05':"../../experiment/PbPb_2p76/charged/ALICE/PbPb_RAA_cent_00_05_eta_0p8.csv",
                   '05-10':"../../experiment/PbPb_2p76/charged/ALICE/PbPb_RAA_cent_05_10_eta_0p8.csv",
                   '10-20':"../../experiment/PbPb_2p76/charged/ALICE/PbPb_RAA_cent_10_20_eta_0p8.csv",
                   '00-20':"../../experiment/PbPb_2p76/charged/ALICE/PbPb_RAA_cent_00_20_eta_0p8.csv"}

raa_fname_cms = {'00-05' : "../../experiment/PbPb_2p76/charged/CMS/Table5_1.csv",
                 '05-10' : "../../experiment/PbPb_2p76/charged/CMS/Table5_2.csv"}

raa_results = {'alice': raa_fname_alice,
               'cms'  : raa_fname_cms,
               'atlas1': raa_fname_atlas1,
               'atlas2': raa_fname_atlas2}

colnames_raa_atlas = {'PT [GeV]':'x',
                      'PT [GeV] LOW' : 'xlow',
                      'PT [GeV] HIGH': 'xhigh',
                      'RAA' : 'y',
                      'stat +' : 'dy_stat+',
                      'stat -' : 'dy_stat-',
                      'sys +' : 'dy_syst+',
                      'sys -' : 'dy_syst-'}
colnames_raa_alice = {'PT [GEV]' : 'x',
                      'PT [GEV] LOW' : 'xlow',
                      'PT [GEV] HIGH': 'xhigh',
                      'RAA' : 'y',
                      'stat +' : 'dy_stat+',
                      'stat -' : 'dy_stat-',
                      'sys +' : 'dy_syst+',
                      'sys -' : 'dy_syst-',
                      "sys,normaliztion uncertainty +" : 'dy_norm+',
                      "sys,normaliztion uncertainty -" : 'dy_norm-'}
colnames_raa_cms = {'PT [GEV]' : 'x',
                    'PT [GEV] LOW' : 'xlow',
                    'PT [GEV] HIGH': 'xhigh',
                    'RAA' : 'y',
                    'stat +' : 'dy_stat+',
                    'stat -' : 'dy_stat-',
                    'sys +' : 'dy_syst+',
                    'sys -' : 'dy_syst-'}

column_names_RAA = {'alice':colnames_raa_alice,
                    'atlas1':colnames_raa_atlas,
                    'atlas2':colnames_raa_atlas,
                    'cms'  :colnames_raa_cms}

available_centralities = {'atlas': ['00-05','05-10','10-20'],
                          'alice': ['00-05','05-10','10-20','00-20'],
                          'cms'  : ['00-05','05-10']} 
from matplotlib.colors import TABLEAU_COLORS as tab

rate_colours = {'LO': '#009e5b', 'NLO':'#002cf1', 'NP':'#ff0000'}
colors = {5:'#e50000', 0:'#d1663c', 3:'#bccc79', 1:'#8fcc9a', 2:'#47669f', 4:'#0000a4'}

rate_colours = {'LO':colors[0], 'NLO':colors[1], 'NP':colors[2]}

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

#colnames_CMS_jet_spectra = {'PT [GEV]' : 'x',
#                'PT [GEV] LOW': 'xlow',
#                'PT [GEV] HIGH': 'xhigh',
#                'E*D3(N)/DP**3 [GEV**-2]' :'y',
#                'stat +' : 'dy_stat+',
#                'stat -' : 'dy_stat-',
#                'sys +' : 'dy_syst+',
#                'sys -': 'dy_syst-'}
colnames_ATLAS_jet_spectra = {'PT [GEV]' : 'x', 
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
colnames_ALICE_jet_spectra = {'PT [GEV/C]' : 'x',
                              'PT [GEV/C] LOW' : 'xlow',
                              'PT [GEV/C] HIGH': 'xhigh',
                              'D2N/DPT/DETA/NEVT [1/(GEV/C)]' : 'y',
                              'stat +' : 'dy_stat+',
                              'stat -' : 'dy_stat-',
                              'sys,uncorr +': 'dy_syst1+',
                              'sys,uncorr -': 'dy_syst1-',
                              'sys,total inelastic cross section uncertainty +' : 'dy_syst2+',
                              'sys,total inelastic cross section uncertainty -' : 'dy_syst2-'}
colnames_CMS_jet_spectra = {'Jet $p_T$ [GEV/c]': 'x',
                            'Jet $p_T$ [GEV/c] LOW': 'xlow',
                            'Jet $p_T$ [GEV/c] HIGH': 'xhigh',
                            '$d^2 \sigma / dp_T d\eta$ [nb/GEV/c]': 'y',
                            'stat +': 'dy_stat+',
                            'stat -': 'dy_stat-',
                            'sys +': 'dy_syst+',
                            'sys -': 'dy_syst-'}
## pp charged hadron spectra:

colnames_ALICE_pp_charged_spectra = {'PT [GEV]' : 'x', 
                                     'PT [GEV] LOW': 'xlow',
                                     'PT [GEV] HIGH': 'xhigh',
                                     '(1/(Nevt))*D2(N)/DETARAP/DPT [C/GEV]': 'y',
                                     'stat +': 'dy_stat+',
                                     'stat -': 'dy_stat-',
                                     'sys +': 'dy_syst+',
                                     'sys -': 'dy_syst-'}
colnames_ATLAS_pp_charged_spectra = {'PT [GeV]' : 'x',
                                     'PT [GeV] LOW' : 'xlow',
                                     'PT [GeV] HIGH': 'xhigh',
                                     '1/(2*PI*PT) D2(SIG)/DPT/DETA [MB/GeV2]' : 'y',
                                     'stat +' : 'dy_stat+',
                                     'stat -' : 'dy_stat-',
                                     'sys +' : 'dy_syst+',
                                     'sys -' : 'dy_syst-'}
colnames_CMS_pp_charged_spectra = {'PT [GEV]' : 'x',
                                   'PT [GEV] LOW': 'xlow',
                                   'PT [GEV] HIGH': 'xhigh',
                                   'E*D3(N)/DP**3 [GEV**-2]' :'y',
                                   'stat +' : 'dy_stat+',
                                   'stat -' : 'dy_stat-',
                                   'sys +' : 'dy_syst+',
                                   'sys -': 'dy_syst-'}
