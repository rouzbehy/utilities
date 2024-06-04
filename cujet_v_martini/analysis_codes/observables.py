"""
    @author: Rouzbeh Yazdi
    @date: July 18 2022
    @description: 
        - Collection of cuts and bins for the varius observables 
        that are used when reading the evolution results.
        - each observable gets its own dictionary
        - the python script calling these dicts is expected to 
        know what's here and to use it properly
"""
from numpy import pi
## pTHat Bins
pT_LHC = [5  ,7  ,9  ,11 ,13 ,15 ,17 ,20 ,\
          25 ,30 ,35 ,40 ,45 ,50 ,55 ,60 ,\
          70 ,80 ,90 ,100,110,120,130,140,\
          150,160,170,180,190,200,210,220,\
          230,240,250,260,270,280,290,300,\
          350,400,450,500,0]

pT_RHIC = [5, 7, 9, 11, 13, 15, 17, 20, 25, 30, 35,\
           40, 45, 50, 55, 60, 65, 70, 75, 80, 85, 90,\
           95, 0]

pT_hat_bins = {0: pT_RHIC, 1: pT_LHC, 2: pT_LHC}

## number of runs per pTHat bin
num_runs_old = 5
num_runs_new = 4
nruns = {'old':num_runs_old, 'new':num_runs_new}

## Charged Hadrons
charged_hadrons_2p76 = {'eta cut': 1, 'spectrum bins':[1,2,3,4,5,6,7,8,9,10,12,14,16,18,20,23,\
                                                     26,30,35,40,45,50,55,60,70,80,90,100,120,\
                                                     140,160,180,200,250,300,350,400,500, 9999]}

charged_hadrons_5p02 = {'eta cut': 1, 'spectrum bins':[1,2,3,4,5,6,7,8,9,10,12,14,16,18,20,23,\
                                                     26,30,35,40,45,50,55,60,70,80,90,100,120,\
                                                     140,160,180,200,250,300,350,400,500, 9999]}

charged_hadrons_0p20 = {'eta cut':0.5, 'spectrum bins':[1,1.5,2.,2.5,3.,3.5,4.,4.5,5.,6,7,8,9,10,\
                                                        12,14,16]}

charged_hadrons = {0:charged_hadrons_0p20, 
                   1:charged_hadrons_2p76,
                   2:charged_hadrons_5p02}

## Photon spectra:
photons_0p20 = {'eta cut':0.35, 'spectrum bins':[4. ,4.5,5. ,5.5,6. ,6.5,7. ,7.5,8.,\
                                                 8.5,9. ,9.5,10.,12.,14.,16.,17.,22]}

photons_2p76 = {'eta cut':0.8, 'spectrum bins':[0.9, 1.1, 1.3, 1.5, 1.7, 1.9, 2.1, 2.3,\
                                                2.5,2.7,3.,3.3,3.7,4.1,4.6,5.4,6.2,7.,\
                                                8.,11.,14.,20.,25.]}

photons_5p02 = {'eta cut':0.8, 'spectrum bins':[0.9, 1.1, 1.3, 1.5, 1.7, 1.9, 2.1, 2.3,\
                                                2.5,2.7,3.,3.3,3.7,4.1,4.6,5.4,6.2,7.,\
                                                8.,11.,14.,20.,25.]}

photon_specs = {0:photons_0p20, 1:photons_2p76, 2:photons_5p02}

## photon v2
photon_v2_0p20 = {'eta cut':0.35, 'spectrum bins':[0.9, 1.1, 1.3, 1.5, 1.7, 1.9, 2.1, 2.3,\
                                                2.5,2.7,3.,3.3,3.7,4.1,4.6,5.4,6.2,7.,\
                                                8.,11.,14.,20.,25.]}

photon_v2_2p76 = {'eta cut':0.8, 'spectrum bins':[0.9, 1.1, 1.3, 1.5, 1.7, 1.9, 2.1, 2.3,\
                                                2.5,2.7,3.,3.3,3.7,4.1,4.6,5.4,6.2,7.,\
                                                8.,11.,14.,20.,25.]}


photon_v2_5p02 = {'eta cut':0.8, 'spectrum bins':[0.9, 1.1, 1.3, 1.5, 1.7, 1.9, 2.1, 2.3,\
                                                2.5,2.7,3.,3.3,3.7,4.1,4.6,5.4,6.2,7.,\
                                                8.,11.,14.,20.,25.]}

photon_v2 = {0:photon_v2_0p20, 1:photon_v2_2p76, 2:photon_v2_5p02}

#oversampling_factors = {'martini_0' : {'conv':10000,'brem':10000}, 
#                        'cujet'     : {'conv':10000,'brem':1000}}
## Jets Spectra:
        
jet_spec_2p76 = {'eta cut':2, 'spectrum bins':[10,12,14,16,18,20,22,25,29,34,40,47,55,64,\
                                                 74,84,97,114,133,153,174,196,220,245,272,\
                                                 300,350,400,450, 500,550,600,700,800,900,1000]}
jet_spec_5p02 = {'eta cut':2, 'spectrum bins':[10,12,14,16,18,20,22,25,29,34,40,47,55,64,\
                                            74,84,97,114,133,153,174,196,220,245,272,\
                                            300,350,400,450, 500,550,600,700,800,900,1000]}
jet_spec_0p20 = {'eta cut':1, 'spectrum bins':[10, 12, 14, 16, 18, 20, 22, 25, 30, 35, 40,\
                                               45, 50, 55, 60, 70, 80, 90, 100]}
jet_spectra = {0:jet_spec_0p20, 1:jet_spec_2p76, 2:jet_spec_5p02}

## for jet fragmentation function
jet_FF_LHC = {'eta cut':2.1, 'pT bins'  :[0.0, 1.0, 1.6, 2.5, 4, 6.3, 10, 16, 25, 40, 63, 100, 9999],
                             'z bins'   :[0.0, 0.01, 0.016, 0.025, 0.04, 0.063, 0.1, 0.16, 0.25, 0.4, 0.63, 1.0],
                             'spec bins':[100,126,158,398]}

jet_FF_RHIC = {'eta cut':1., 'pT bins'  :[0.0, 1.0, 1.6, 2.5, 4, 6.3, 10, 16, 25, 40, 63, 100],
                             'z bins'   :[0.0, 0.01, 0.016, 0.025, 0.04, 0.063, 0.1, 0.16, 0.25, 0.4, 0.63, 1.0],
                             'spec bins':[10,30,60,90]}
                    
## for jet shape function
jet_shape_LHC = {'eta min' : 0.3, 'eta max':2, 'pT min':100, 'pT max':9999, 'r bins':[0, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35]}
jet_shape_RHIC = {'eta min' : -1.0, 'eta max':1.0, 'pT min':10, 'pT max': 100, 'r bins':[0, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35]}

## for dijet observables
PI = 3.141592653589793
dijet_hard_prop = {'eta cut': 2.1, 'phi min':7*pi/8.,'xj bins':[0.32, 0.35, 0.40, 0.45, 0.50, 0.56, 0.63, 0.71, 0.79, 0.89, 1.00]}
dijet_hard_pT = [100,112,126,141,158,178,200,224,251,282,316,398,562,627,700]

dijet_minijet_prop = {'eta cut': 1, 'phi min':7*pi/8., 'xj bins':[0.32, 0.35, 0.40, 0.45, 0.50, 0.56, 0.63, 0.71, 0.79, 0.89, 1.00]}
dijet_minijet_pT = [10, 14, 18, 22, 26, 30, 34, 40, 46, 52, 58, 65, 73, 82, 92, 100]
