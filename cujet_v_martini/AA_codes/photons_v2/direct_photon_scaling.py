import numpy as np
import pandas as pd
from scipy.optimize import curve_fit
from scipy.interpolate import InterpolatedUnivariateSpline
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import util
import jetDicts as ddicts
from dictionaries import multiplicity_AuAu_200, multiplicity_PbPb_2p76, multiplicity_PbPb_5p02
from dictionaries import oversampling_factors as osfact_2p76
from COLORS import module_colors, colors, channel_colors
import photonfuncs as pfuncs

plt.rcParams.update(util.my_rcParams)

## The scaling function to fit to
def func(x, logA, alpha):
    return logA + x*alpha

## Read in the charged hadron multiplicity file
ch_mult = pd.read_csv('ch_mult.csv')
ch_mult.sort_values(by=['nch'], inplace=True)
print(ch_mult)
cross_sections = {'AuAu200':42, 'PbPb2760':62, 'PbPb5020':68}
centralities = {'AuAu200':['00-10','10-20'],
                'PbPb2760':['00-05','05-10','10-20','20-30','30-40','40-50'],
                'PbPb5020':['00-10','10-20','30-50']}
kfactors = {'AuAu200':1.13,'PbPb2760':0.71,'PbPb5020':0.9}
multiplicities = {'AuAu200':multiplicity_AuAu_200, 
                  'PbPb2760':multiplicity_PbPb_2p76,
                  'PbPb5020':multiplicity_PbPb_5p02}
dumb_name = {'AuAu200':'AuAu_200', 'PbPb2760':'PbPb_2760', 'PbPb5020':'PbPb_5020'}

ROOT_DIRECTORY   = "/Users/rmyazdi/Documents/research/jetscape_project/"
THERMAL_ROOT_DIR = "/Users/rmyazdi/Documents/research/hydro_plots/VISHNU/HydroProfiles/"
KOMPOST_ROOT_DIR = ROOT_DIRECTORY + 'v2/other_data/JF_MultiMessenger/'
PROMPT_ROOT_DIR  = ROOT_DIRECTORY+"v2/jetscape_data/prompt_photons/"
JMEDIUM_ROOT_DIR = ROOT_DIRECTORY + 'v2/jetscape_data/'
## read in the data, construct three sets:
martini, cujet, prompt, kompost, thermal = {}, {}, {}, {}, {}
pT_min_values = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18]
pT_max = 1000
two_pi = 2 * np.pi
for pT_min in pT_min_values:
    martini[pT_min] = np.zeros_like(ch_mult['nch']) 
    cujet[pT_min]   = np.zeros_like(ch_mult['nch']) 
    prompt[pT_min]  = np.zeros_like(ch_mult['nch']) 
    thermal[pT_min] = np.zeros_like(ch_mult['nch']) 
    kompost[pT_min] = np.zeros_like(ch_mult['nch']) 
    for index, row in ch_mult.iterrows():
        nuclei, energy, centrality = row['nuclei'], row['energy'], row['centrality']
        tag1 = f'{nuclei}{energy}'
        xsec = cross_sections[tag1]
        multiplicity = multiplicities[tag1]
        kfactor = kfactors[tag1]

        ## channel by channel:
        loc = PROMPT_ROOT_DIR+dumb_name[tag1]+'/'
        tmp, _ = pfuncs.get_prompt_photons(loc, multiplicity[centrality], None, xsec, kfactor=1)
        tmp = tmp[tmp['pT'].between(pT_min, pT_max)]
        val = sum([two_pi * dx * x * y for (x,dx,y) in zip(tmp['pT'],tmp['dpT'],tmp['N'])])
        prompt[pT_min][index]=val
        #print(f"prompt: {pT_min}, {val}")

        tmp = pfuncs.get_JF_photons(tag1, centrality, 'preEq', KOMPOST_ROOT_DIR)
        tmp = tmp[tmp['pT'].between(pT_min, pT_max)]
        val = sum([two_pi * 0.2 * x * y for (x,y) in zip(tmp['pT'],tmp['N'])]) 
        kompost[pT_min][index]=val
        #print(f"kompost: {pT_min}, {val}")

        tmp = pfuncs.get_JF_photons(tag1, centrality, 'thermal', KOMPOST_ROOT_DIR)#pfuncs.get_thermal_photons(tag1, centrality, THERMAL_ROOT_DIR)
        tmp = tmp[tmp['pT'].between(pT_min, pT_max)]
        val = sum([two_pi * 0.2 * x * y for (x,y) in zip(tmp['pT'],tmp['N'])])
        #f = InterpolatedUnivariateSpline(y=two_pi*tmp['pT']*tmp['N'], x=tmp['pT'])
        #val = f.integral(pT_min, pT_max)
        thermal[pT_min][index]=val
        #print(f"thermal: {pT_min}, {val}")

        tmp = pfuncs.get_jet_medium_photons('martini', centrality, xsec, energy, tag1, multiplicity, JMEDIUM_ROOT_DIR)
        tmp = tmp[tmp['pT'].between(pT_min, pT_max)]
        val = sum([two_pi*x*dx*y for (x,y, dx) in  zip(tmp['pT'], tmp['total'], tmp['dpT'])])
        martini[pT_min][index]=val
        #print(f"martini: {pT_min}, {val}")

        tmp = pfuncs.get_jet_medium_photons('cujet', centrality, xsec, energy, tag1, multiplicity, JMEDIUM_ROOT_DIR)
        tmp = tmp[tmp['pT'].between(pT_min, pT_max)]
        val = sum([two_pi*x*dx*y for (x,y, dx) in  zip(tmp['pT'], tmp['total'], tmp['dpT'])])
        cujet[pT_min][index]=val
        #print(f"cujet: {pT_min}, {val}")

## Construct the totals

totals = {'martini':{}, 'cujet':{}, 'no jet-med.':{}}
template = ','.join(["{:0.1e}" for item in range(len(ch_mult['nch']))])
parameters = {'martini':{},'cujet':{},'no jet-med.':{}}
alphas, Avalues = [], []
for pT_min in pT_min_values:
    no_jet_med = prompt[pT_min] + thermal[pT_min] + kompost[pT_min]
    martini_spec = no_jet_med + martini[pT_min]
    cujet_spec = no_jet_med + cujet[pT_min]
    
    # print(f"{pT_min}: ")
    # print("\t no jet med: " + template.format(*no_jet_med))
    # print("\t martini   : " + template.format(*martini_spec)) 
    # print("\t cujet     : " + template.format(*cujet_spec))
    # print("---------------------------------")
    totals['cujet']      [pT_min] = cujet_spec
    totals['martini']    [pT_min] = martini_spec
    totals['no jet-med.'][pT_min] = no_jet_med

    params, pcov = curve_fit(func, np.log(ch_mult['nch']), np.log(no_jet_med), p0=(1e-3, 1.25))
    parameters['no jet-med.'][pT_min] = params 

    params, pcov = curve_fit(func, np.log(ch_mult['nch']), np.log(martini_spec), p0=[1e-3,1.25])
    parameters['martini'][pT_min] = params

    params, pcov = curve_fit(func, np.log(ch_mult['nch']), np.log(cujet_spec), p0=[1e-3,1.25])
    parameters['cujet'][pT_min] = params

vfunc = np.vectorize(func)
## Now plot it!
fig1, ax1 = plt.subplots(1,1) ## plot of photons vs charged hadrons
fig2, ax2 = plt.subplots(2,1) ## plot of the scaling parameter alpha
for key in totals:
    for pT in [1, 2, 16]:
        params = parameters[key][pT]
        print(f"{key}, {pT}, {np.exp(params[0]):0.5e}, {params[1]:0.3f}")
        ax1.scatter(ch_mult['nch'], totals[key][pT])
        ax1.plot(ch_mult['nch'], np.exp(vfunc(np.log(ch_mult['nch']), *params)))

for key in parameters:
    y1 = [parameters[key][pT][1] for pT in pT_min_values]
    y2 = [np.exp(parameters[key][pT][0]) for pT in pT_min_values] 
    ax2[0].plot(pT_min_values, y1, label=key)
    ax2[1].plot(pT_min_values, y2, label=key)

ax2[0].legend(loc='best')
ax2[1].set_yscale('log')
#ax1.set_yscale('log')
#ax1.set_xscale('log')
ax1.set_ylim(top=54,bottom=1e-12)
plt.show()