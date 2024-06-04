#!/usr/bin/env python3
"""
    Plot the computed spectra for jets 
"""
## import scipy and matplotlib modules:
import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
## my custom modules
import util 
import jetDicts as ddicts
import dictionaries as my_dicts
#mpl.use('Qt5Agg')
my_rcParams={
    "text.usetex": True,
    "font.family": "Georgia",
    "font.size": 25,
    "lines.linewidth": 2,
    "xtick.direction": "in",
    "ytick.direction": "in",
    "xtick.minor.visible": True,
    "ytick.minor.visible": True,
    "xtick.major.size" : 8,
    "ytick.major.size" : 8,
    "xtick.minor.size" : 4,
    "ytick.minor.size" : 4,
    "axes.spines.right": False,
    "axes.spines.top" : False,
    "legend.frameon":False
}
plt.rcParams.update(my_rcParams)

## read in jetscape calculation
minpT, maxpT = 0.1, 120
zmin, zmax = 0.005,2

def get_FF_martini(fname, fftype):
    tmp = pd.read_csv(fname, comment='#')
    njets, nevts = 0 , 0
    if fftype == 'pT':
        tmp = tmp[tmp['pTmax'] < maxpT]
        tmp = tmp[tmp['pTmin'] > minpT]
        dx = tmp['pTmax'] - tmp['pTmin']
        x = 0.5*(tmp['pTmax'] + tmp['pTmin'])
    if fftype == 'z':
        tmp = tmp[tmp['zmax'] < zmax]
        tmp = tmp[tmp['zmin'] > zmin]
        dx = tmp['zmax'] - tmp['zmin']
        x = 0.5*(tmp['zmax'] + tmp['zmin'])
    with open(fname,'r') as f:
        line = f.readline()
        line = line.split(' ')
        njets = float(line[1])
        nevts = float(line[-1])
    scale = njets
    y  = tmp['N'] /(dx*scale)
    dy = tmp['dN']/(dx*scale)   
    y  = tmp['N'] /(dx*scale)
    dy = tmp['dN']/(dx*scale)
    result = pd.DataFrame({'x':x,'dx':dx, 'y':y, 'dy':dy})
    return result

def get_FF_jetscape(directory, system='PbPb2760'):
        if system == 'PbPb2760':
            fname_pT   = directory + f'PbPb2760_00-05_LHC_FF_pT_jet_rad_0.4_2.10.csv'
            fname_z    = directory + f'PbPb2760_00-05_LHC_FF_z_jet_rad_0.4_2.10.csv'
            fname_spec = directory + f'PbPb2760_00-05_LHC_jet_spec_FF_jet_rad_0.4_2.10.csv'
        else:
            fname_pT   = directory + f'pp_2760_LHC_FF_pT_jet_rad_0.4_2.10.csv'
            fname_z    = directory + f'pp_2760_LHC_FF_z_jet_rad_0.4_2.10.csv'
            fname_spec = directory + f'pp_2760_LHC_jet_spec_FF_jet_rad_0.4_2.10.csv'

        tmp_pT   = pd.read_csv(fname_pT, comment='#')
        tmp_z    = pd.read_csv(fname_z,  comment='#')
        tmp_spec = pd.read_csv(fname_spec, comment='#')
        njets = sum(tmp_spec['Ncut'].tolist())
        tmp_pT = tmp_pT[tmp_pT['pTmax'] < maxpT] 
        tmp_pT = tmp_pT[tmp_pT['pTmin'] > minpT]
        dpT = tmp_pT['pTmax'] - tmp_pT['pTmin'] 
        pT  = 0.5*(tmp_pT['pTmax'] + tmp_pT['pTmin'])
        y_pT = tmp_pT['Ncut']/(njets * dpT)
        dy_pT = tmp_pT['dNcut']/(njets * dpT)
        tmp_z = tmp_z[tmp_z['zmax'] < zmax] 
        tmp_z = tmp_z[tmp_z['zmin'] > zmin]
        dz = tmp_z['zmax'] - tmp_z['zmin'] 
        z  = 0.5*(tmp_z['zmax'] + tmp_z['zmin'])
        y_z  = tmp_z['Ncut']/(njets * dz)
        dy_z = tmp_z['dNcut']/(njets * dz)
        result_pT = pd.DataFrame({'x':pT,'dx':dpT,'y':y_pT,'dy':dy_pT})
        result_z  = pd.DataFrame({'x':z, 'dx':dz, 'y':y_z, 'dy':dy_z})
        return result_pT, result_z

from COLORS import colors, shower_effects

if __name__=='__main__':
    ## Read in MARTINI Stand-Alone work:
    martini_pp_DPT = get_FF_martini('../martini_results/pp/jet_FF_pT.csv','pT')
    martini_pp_DZ  = get_FF_martini('../martini_results/pp/jet_FF_z.csv','z')
    martini_AA_DPT = get_FF_martini('../martini_results/final_PbPb_2p76/rset_1/cent_0_5/jet_FF_pT.csv', 'pT')
    martini_AA_DZ  = get_FF_martini('../martini_results/final_PbPb_2p76/rset_1/cent_0_5/jet_FF_z.csv', 'z') 

    ## Read in JETSCAPE MARTINI work
    pp_dir = '../../../../jetscape_project/v2/jetscape_data/max_time/maxT_200_highstat/'
    jetscape_pp_DPT, jetscape_pp_DZ = get_FF_jetscape(pp_dir, system='pp_2760')
    old_kappa_dir = '../../../../jetscape_project/v2/jetscape_data/sqrt_s_2760/martini/PbPb_2760/'
    jetscape_AA_DPT_old, jetscape_AA_DZ_old = get_FF_jetscape(old_kappa_dir)
    new_kappa_dir = '../../../../jetscape_project/v2/jetscape_data/sqrt_s_2760/martini_new_kappas/'
    jetscape_AA_DPT_new, jetscape_AA_DZ_new = get_FF_jetscape(new_kappa_dir)
    new_kappa_dir_2 = '../../../../jetscape_project/v2/jetscape_data/sqrt_s_2760/martini_new_kap_matter_vac/'
    jetscape_AA_DPT_new_vac, jetscape_AA_DZ_new_vac = get_FF_jetscape(new_kappa_dir_2)


    exp_loc = '../../../exp_data/jet_fragmentation/ATLAS/'
    ## read in ATLAS results
    data_pT = pd.read_csv(exp_loc+'HEPData-ins1511869-v1-csv/Table9.csv',comment='#').rename(columns=ddicts.colnames_ATLAS_FF_DPT)
    data_z = pd.read_csv(exp_loc+"HEPData-ins1511869-v1-csv/Table25.csv",comment='#').rename(columns=ddicts.colnames_ATLAS_FF_DZ)
    ## make the plots:
    fig, axes = plt.subplots(1, 2, figsize=(16,9), sharex='col')

    util.plot_expr_data_on_axis(axes[0], data_pT, 's')

    ## MARTINI only curve
    pT = martini_pp_DPT['x']
    yp1, dyp1 = martini_pp_DPT['y'], martini_pp_DPT['dy']
    ya1, dya1 = martini_AA_DPT['y'], martini_AA_DPT['dy']
    raa1 = ya1/yp1
    draa1 = raa1*np.sqrt(dyp1*dyp1/(yp1*yp1) + dya1*dya1/(ya1*ya1))
    axes[0].plot(pT, raa1, color=colors[0])
    axes[0].fill_between(pT, raa1-draa1, raa1+draa1, color=colors[0], alpha=0.2)


    pT = jetscape_pp_DPT['x']
    yp1, dyp1 = jetscape_pp_DPT['y'], jetscape_pp_DPT['dy']
    ya1, dya1 = jetscape_AA_DPT_old['y'], jetscape_AA_DPT_old['dy']
    raa1 = ya1/yp1
    draa1 = raa1*np.sqrt(dyp1*dyp1/(yp1*yp1) + dya1*dya1/(ya1*ya1))
    axes[0].plot(pT, raa1, color=colors[1])
    axes[0].fill_between(pT, raa1-draa1, raa1+draa1, color=colors[1], alpha=0.2)

    ya1, dya1 = jetscape_AA_DPT_new['y'], jetscape_AA_DPT_new['dy']
    raa1 = ya1/yp1
    draa1 = raa1*np.sqrt(dyp1*dyp1/(yp1*yp1) + dya1*dya1/(ya1*ya1))
    axes[0].plot(pT, raa1, color=colors[2])
    axes[0].fill_between(pT, raa1-draa1, raa1+draa1, color=colors[2], alpha=0.2)

    ya1, dya1 = jetscape_AA_DPT_new_vac['y'], jetscape_AA_DPT_new_vac['dy']
    raa1 = ya1/yp1
    draa1 = raa1*np.sqrt(dyp1*dyp1/(yp1*yp1) + dya1*dya1/(ya1*ya1))
    axes[0].plot(pT, raa1, color=colors[3])
    axes[0].fill_between(pT, raa1-draa1, raa1+draa1, color=colors[3], alpha=0.2)

    ## MARTINI only curve
    util.plot_expr_data_on_axis(axes[1], data_z, 's')
    pT = martini_pp_DZ['x']
    yp1, dyp1 = martini_pp_DZ['y'], martini_pp_DZ['dy']
    ya1, dya1 = martini_AA_DZ['y'], martini_AA_DZ['dy']
    raa1 = ya1/yp1
    draa1 = raa1*np.sqrt(dyp1*dyp1/(yp1*yp1) + dya1*dya1/(ya1*ya1))
    axes[1].plot(pT, raa1, color=colors[0])
    axes[1].fill_between(pT, raa1-draa1, raa1+draa1, color=colors[0], alpha=0.2)

    pT = jetscape_pp_DZ['x']
    yp1, dyp1 = jetscape_pp_DZ['y'], jetscape_pp_DZ['dy']
    ya1, dya1 = jetscape_AA_DZ_old['y'], jetscape_AA_DZ_old['dy']
    raa1 = ya1/yp1
    draa1 = raa1*np.sqrt(dyp1*dyp1/(yp1*yp1) + dya1*dya1/(ya1*ya1))
    axes[1].plot(pT, raa1, color=colors[1])
    axes[1].fill_between(pT, raa1-draa1, raa1+draa1, color=colors[1], alpha=0.2)

    ya1, dya1 = jetscape_AA_DZ_new['y'], jetscape_AA_DZ_new['dy']
    raa1 = ya1/yp1
    draa1 = raa1*np.sqrt(dyp1*dyp1/(yp1*yp1) + dya1*dya1/(ya1*ya1))
    axes[1].plot(pT, raa1, color=colors[2])
    axes[1].fill_between(pT, raa1-draa1, raa1+draa1, color=colors[2], alpha=0.2)

    ya1, dya1 = jetscape_AA_DZ_new_vac['y'], jetscape_AA_DZ_new_vac['dy']
    raa1 = ya1/yp1
    draa1 = raa1*np.sqrt(dyp1*dyp1/(yp1*yp1) + dya1*dya1/(ya1*ya1))
    axes[1].plot(pT, raa1, color=colors[3])
    axes[1].fill_between(pT, raa1-draa1, raa1+draa1, color=colors[3], alpha=0.2)

    labels = [Line2D([],[],color=c,label=l) for l, c in shower_effects.items()]
    #labels.append()
    #axes[0].legend(loc='upper left', handles=labels, fontsize=18)
    artist= axes[0].legend(handles=labels, fontsize=18, bbox_to_anchor=(0.1,1.1), loc='upper left' )
    axes[0].add_artist(artist)
    axes[0].set_xscale('log')
    axes[0].set_ylabel(r'$R_{D(p_T)}$')
    axes[1].set_ylabel(r'$R_{D(z)}$')
    axes[0].set_xlabel(r'$p^{h^{\pm}}_T$ (GeV)')
    axes[1].set_xlabel(r'$z$')

    expt_label = [Line2D([],[],color='black', marker='s', label='ATLAS (2017) $0$-$10\%$')]
    axes[1].legend(loc='upper left', handles=expt_label)
    axes[1].text(0.1,0.8,s=r'Pb-Pb, $\sqrt{s}=2.76$ ATeV'         , transform=axes[1].transAxes, fontsize=18)
    axes[1].text(0.1,0.75,s=r'$100< p^{\mathrm{jet}}_T < 398$ GeV', transform=axes[1].transAxes, fontsize=18)
    axes[1].text(0.1,0.7,s=r'$0$-$5\%$'                          , transform=axes[1].transAxes, fontsize=18)
    plt.show()
