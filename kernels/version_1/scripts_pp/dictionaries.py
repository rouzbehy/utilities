from numpy import sqrt
multiplicity = {'00-05':1615,
                '05-10':1268,
                '10-15':1002,
                '15-20':790.1,
                '20-30':548.5,
                '30-40':320.9,
                '40-50':174.85}
#'00-20':(1615+1268+896.05)/4,
#'20-40':(618.1 +478.9+366.3+275.5)/4}
# multiplicity = {'00-05':1615+135,
#                 '00-20':(1615+1268+896.05)/4,
#                 '20-30': (618.1+66 + 478.9+59.1)/2,
#                 '30-40': (366.3+52.7 + 375.5+46.7)/2}
#multiplicity_err = {'00-05':135,
#                    '20-30':sqrt(66**2 + 59.1**2)/2,
#                    '30-40':sqrt(52.7**2 + 46.7**2)/2,
#                    '20-40':sqrt(66**2 + 59.1**2 + 52.7**2 + 46.7**2)/2}
exag_factor = {'cujet'    : 500000,
               'martini_0': 1000000,
               'martini_1': 1000000}

cent_colours = {'00-05' : "#0077BB",
                '20-40' : "#CC3311",
                '20-30' : "#33BBEE",
                '30-40' : "#009988",
                }
eloss_color = {
    'cujet' : "#CC3311",
    'martini_0':"#0077BB",
    'martini_1':"tab:olive",
    'old Martini':"#009988"
}
centrality_lstyle = {'00-05': 'solid',
                     '30-40': 'dashed'}

channel_color = {'conv-cujet'    :"#CC3311",
                 'conv-martini':"#0077BB",
                 'prompt'        :"#33BBEE",
                 'thermal'       :"#009988",
                 'brem-martini':"#EE3377"}


nice_cent_label = {'00-05': '0-5\%', 
                   '20-30': '20-30\%',
                   '30-40': '30-40\%',
                   '20-40': '20-40\%'}
JF_cent_labels = {"C0-5":"00-05","C20-30":"20-30","C30-40":"30-40"}

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
colnames_raa_alice = {
    "PT [GEV]":"x",
    "PT [GEV] LOW":"xlow",
    "PT [GEV] HIGH":"xhigh",
    "RAA":"y",
    "stat +":"dy_stat+",
    "stat -":"dy_stat-",
    "sys +":"dy_syst+",
    "sys -":"dy_syst-",
    "sys,normaliztion uncertainty +":"dy_norm+",
    "sys,normaliztion uncertainty -":"dy_norm-"
}
## pp charged hadron spectra:

colnames_ALICE_pp_chHad_spec = {'PT [GEV]' : 'x', 
                                     'PT [GEV] LOW': 'xlow',
                                     'PT [GEV] HIGH': 'xhigh',
                                     '(1/(Nevt))*D2(N)/DETARAP/DPT [C/GEV]': 'y',
                                     'stat +': 'dy_stat+',
                                     'stat -': 'dy_stat-',
                                     'sys +': 'dy_syst+',
                                     'sys -': 'dy_syst-'}
colnames_ATLAS_pp_chHad_spec = {'PT [GeV]' : 'x',
                                     'PT [GeV] LOW' : 'xlow',
                                     'PT [GeV] HIGH': 'xhigh',
                                     '1/(2*PI*PT) D2(SIG)/DPT/DETA [MB/GeV2]' : 'y',
                                     'stat +' : 'dy_stat+',
                                     'stat -' : 'dy_stat-',
                                     'sys +' : 'dy_syst+',
                                     'sys -' : 'dy_syst-'}
colnames_CMS_pp_chHad_spec = {'PT [GEV]' : 'x',
                                   'PT [GEV] LOW': 'xlow',
                                   'PT [GEV] HIGH': 'xhigh',
                                   'E*D3(N)/DP**3 [GEV**-2]' :'y',
                                   'stat +' : 'dy_stat+',
                                   'stat -' : 'dy_stat-',
                                   'sys +' : 'dy_syst+',
                                   'sys -': 'dy_syst-'}


## RHIC Crap
colnames_STAR_pp_ch_spec = {'$p_T [GeV/c]$':'x',
                            '$p_T [GeV/c]$ LOW':'xlow',
                            '$p_T [GeV/c]$ HIGH':'xhigh',
                            "$1/(2\pi p_T) d^2N/dpTd\eta, |\eta|<0.5 [(GeV/c)^-2]$":'y',
                            'sys + stat +':'dy_+',
                            'sys + stat -':'dy_-'}