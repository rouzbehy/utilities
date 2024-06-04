import os
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.interpolate import InterpolatedUnivariateSpline as interp1d
import dictionaries as hDicts
import jetDicts as jDicts
import util

plt.rcParams.update(util.my_rcParams)

def calc_Xi_sq(data, f_th, df_th,  include_theory_err=False):
    tx = data['x']
    ty = f_th(tx)
    dty = df_th(tx)
    chisq = (ty - data['y'])**2 #array([(v1-v2)**2 for (v1,v2) in zip(ty,data['y'])])
    errs  = (data['dy_stat+']**2+ data['dy_syst+']**2)#array([v1*v1 + v2*v2 for v1, v2 in zip(data['dy_stat+'], data['dy_syst+'])])
    if include_theory_err:
        errs = errs + dty**2 #array([v1 + e*e for (v1,e) in zip(errs,dty)])
    chisq /= errs
    return sum(chisq), len(data['y'])

def read_in_experiments():
    ## Experimental Data Read-In
    alice_hadron_name  = "../../../exp_data/charged/ALICE/Table_16.csv"
    cms_hads_name      = "../../../exp_data/charged/CMS/Table5_1.csv"
    cms_jets_name      = "../../../exp_data/jets/CMS/Table_JETRAA_0p{i}.csv"
    cms_jet_shape_name = '../../../exp_data/jet_shape/CMS_PbPb_2760_Jet_ShapeRatio_00-10_pT=100-INF_R=0p3_eta=0p3-2p0.dat'
    atlas_hads_name_2  = "../../../exp_data/charged/ATLAS/Table33.csv"
    atlas_hads_name_1  = '../../../exp_data/charged/ATLAS/Table49.csv'
    ALICE_HADS         = pd.read_csv(alice_hadron_name, 
                            comment='#').rename(columns=hDicts.colnames_raa_alice)
    ATLAS_HADS_2       = pd.read_csv(atlas_hads_name_2, 
                            comment='#').rename(columns=hDicts.colnames_raa_atlas)
    ATLAS_HADS_1       = pd.read_csv(atlas_hads_name_1, 
                            comment='#').rename(columns=hDicts.colnames_raa_atlas)
    CMS_HADS           = pd.read_csv(cms_hads_name, 
                            comment='#').rename(columns=hDicts.colnames_raa_cms)
    CMS_JETS           = {f'0p{i}': pd.read_csv(cms_jets_name.format(i=i), 
                            comment='#').rename(columns=jDicts.colnames_CMS_RAA) for i in [2,3,4]}
    tmp                = np.loadtxt(cms_jet_shape_name,comments='#',unpack=True, delimiter='\t')
    CMS_JET_SHAPE      = pd.DataFrame({'x':tmp[0],'y':tmp[1],'dx':tmp[2],
                                    'dy_stat+':tmp[3], 'dy_stat-':tmp[3], 
                                    'dy_syst-':tmp[1]-tmp[4], 'dy_syst+':tmp[5]-tmp[1]})
    return {"ALICE_HADS":ALICE_HADS, "ATLAS_HADS_1":ATLAS_HADS_1, 
            "ATLAS_HADS_2":ATLAS_HADS_2, "CMS_HADS":CMS_HADS, 
            "CMS_JETS":CMS_JETS, "CMS_JET_SHAPE":CMS_JET_SHAPE}

if __name__=='__main__':    
    root_chisq = '../../calcs/chisq_values/'
    root_calcs = '../../calcs/RAA_for_fit/'
    if not os.path.exists(root_chisq):
        os.system(f"mkdir -p {root_chisq}")
    fname_hads       = root_calcs+ 'chgd_raa_rset_{r}_kset_{k}.csv'
    fname_jets       = root_calcs+ 'jet_raa_rset_{r}_kset_{k}.csv'
    fname_jet_shapes = root_calcs+ 'jet_shape_ratio_rset_{r}_kset_{k}.csv'

    rate_sets = util.rate_set_names
    klist     = util.kappas

    pt_chgd_low = util.pt_chgd_low
    pt_chgd_high = util.pt_chgd_high
    pt_jet_low = util.pt_jet_low
    pt_jet_high = util.pt_jet_high
    jet_cone_radii = util.jet_cone_radii
    raw_experiments = read_in_experiments()
    experiments = {}
    for experiment in raw_experiments:
        if 'SHAPE' in experiment:
            experiments[experiment] = raw_experiments[experiment]
            continue
        spec = raw_experiments[experiment]
        if 'HADS' in experiment:
            xlow, xhigh = (pt_chgd_low, pt_chgd_high) 
            experiments[experiment] = spec[spec['x'].between(xlow, xhigh)]
        else:
            xlow, xhigh = (pt_jet_low, pt_jet_high)
            if experiment not in experiments:
                experiments[experiment] = {}
            for r in jet_cone_radii:
                experiments[experiment][r] = spec[r][spec[r]['x'].between(xlow, xhigh)]
        

    calc_hadron    = {rset: {kset: pd.read_csv(fname_hads.format(r=rset,k=kset), comment='#') for kset in np.arange(0,50)} for rset in rate_sets}
    calc_jets      = {rset: {kset: pd.read_csv(fname_jets.format(r=rset,k=kset), comment='#') for kset in np.arange(0,50)} for rset in rate_sets}
    calc_jet_shape = {rset: {kset: pd.read_csv(fname_jet_shapes.format(r=rset,k=kset), comment='#') for kset in np.arange(0,50)} for rset in rate_sets}

    cone_radii = ['0p2','0p3','0p4']
    hadrons = {'ALICE': experiments['ALICE_HADS'], 'CMS':experiments['CMS_HADS'],
                'ATLAS1':experiments['ATLAS_HADS_1'], 'ATLAS2':experiments['ATLAS_HADS_2']}

    fit_data_sets = {'jets': True, 'charged hadrons': True, 'jet shape': False}
    for r in rate_sets:
        print(util.bcolors.OKBLUE, f'Rate Set: {r}', util.bcolors.ENDC)
        results = {'k':[],'kappa_r':[],'kappa_e':[],'chi_sq/ndf':[]}
        for k in np.arange(0,50):
            results['k'].append(k)
            kapr, kape = klist[k]
            results['kappa_r'].append(kapr)
            results['kappa_e'].append(kape)

            ##print('\t',util.bcolors.OKGREEN, f'K set: {k}', util.bcolors.ENDC)
            
            npoints = 0            
            ## CHARGED HADRON SPECTRA
            chisq = 0
            if fit_data_sets['charged hadrons']:
                spec = calc_hadron[r][k]
                pT, raa, draa = spec['pT'], spec['raa'], spec['draa']
                f_hd  = interp1d(pT, raa , k=1)
                df_hd = interp1d(pT, draa, k=1)
                for experiment_name in hadrons:
                    if experiment_name not in results:
                        results[experiment_name] = []
                    dat = hadrons[experiment_name]
                    tmp, n1 = calc_Xi_sq(dat, f_hd, df_hd, True)
                    results[experiment_name].append(tmp)
                    chisq += tmp
                    npoints += n1

            ### JET SPECTRA
            if fit_data_sets['jets']:
                spec    = calc_jets[r][k]
                data = experiments['CMS_JETS']
                for R in cone_radii:
                    if f'CMS_JETS_{R}' not in results:
                        results[f'CMS_JETS_{R}'] = []
                    pT, raa, draa = [spec[label] for label in ['pT',f'r{R}', f'dr{R}']] 
                    f_jet  = interp1d(pT, raa,  k=1)
                    df_jet = interp1d(pT, draa, k=1)
                    tmp, n2 = calc_Xi_sq(data[R], f_jet, df_jet, True)
                    results[f'CMS_JETS_{R}'].append(tmp)
                    chisq += tmp
                    npoints += n2

            ## JET SHAPE (Calculation is for 0-5% , data is for 0-10%, be aware.)  
            if fit_data_sets['jet shape']:
                spec  = calc_jet_shape[r][k]
                data = experiments['CMS_JET_SHAPE']
                x, rho, drho = spec['r'], spec['Rrho'], spec['dRrho']
                if 'CMS_JET_SHAPE' not in results:
                    results['CMS_JET_SHAPE'] = []
                f_shape  = interp1d(x, rho,  k=1)
                df_shape = interp1d(x, drho, k=1)
                tmp, n3 = calc_Xi_sq(data, f_shape, df_shape, include_theory_err=True)
                results['CMS_JET_SHAPE'].append(tmp)
                npoints += n3
                chisq += tmp

            # print(chisq, npoints)
            chisq /= npoints
            results['chi_sq/ndf'].append(chisq)
        results = pd.DataFrame(results)
        print(results.to_markdown())
        results.to_csv(root_chisq+f"rset_{r}.csv", index=False)