import pandas as pd
import numpy as np
from scipy.interpolate import interp1d
from matplotlib.collections import PatchCollection
from matplotlib.patches import Rectangle
import util
from dictionaries import multiplicity_PbPb_2p76 as multiplicity
from dictionaries import oversampling_factors

def get_experimental_data():
    ## read in experimental data
    exp_loc = "../../../expt/PbPb_2p76/photons/"
    centralities = {'00-20':exp_loc+'HEPData-ins1394677-v1-Table_1.csv',
                    '20-40':exp_loc+'HEPData-ins1394677-v1-Table_2.csv'}

    expdata = {}
    pT_lower_lim = {'00-20':0,'20-40':0}
    pT_upper_lim = {'00-20':0,'20-40':0}
    for cent in centralities:
        tmp = pd.read_csv(centralities[cent],comment='#')
        tmp['dx'] = tmp['xhigh'] - tmp['xlow']
        #min_pT = min(tmp['x'])
        max_pT = max(tmp['x'])
        pT_lower_lim[cent] = 2
        pT_upper_lim[cent] = max_pT
        tmp = tmp[tmp['x'].between(pT_lower_lim[cent],pT_upper_lim[cent])]
        expdata[cent] = tmp
    return expdata, pT_lower_lim, pT_upper_lim

def get_jetscape_data(expdata, lower_bound, upper_bound, xsec=1):
    ## read jetscape:
    elosses = ['martini','cujet']
    jetscape_cents = ['00-05','05-10','10-20','20-30','30-40']
    jetscape = {}
    fname_tmpl = '../../jetscape_data/sqrt_s_2760/{eloss}/PbPb_2760/PbPb2760_{cent}_photon_spec_0.80.csv'
    for eloss in elosses:
        jetscape[eloss] = {}
        for cent in jetscape_cents:
            fname = fname_tmpl.format(eloss=eloss,cent=cent)
            factor = oversampling_factors[f'{eloss}-{cent}']
            Nbin = multiplicity[cent]
            tmp = pd.read_csv(fname,comment='#')
            tmp['pT'] = 0.5* (tmp['ptmin'] + tmp['ptmax'])
            tmp['dpT'] = tmp['ptmax'] - tmp['ptmin']
            tmp = tmp[tmp['pT'].between(lower_bound['00-20'],upper_bound['00-20'])]
            tmp['conv']  *= Nbin/(factor*(tmp['pT']*tmp['dpT']*xsec*2*np.pi*2*0.8))
            tmp['dconv'] *= np.sqrt(Nbin)/(np.sqrt(factor)*(tmp['pT']*tmp['dpT']*xsec*2*np.pi*2*0.8))
            tmp['brem']  *= Nbin/(factor*(tmp['pT']*tmp['dpT']*xsec*2*np.pi*2*0.8))
            tmp['dbrem'] *= np.sqrt(Nbin)/(np.sqrt(factor)*(tmp['pT']*tmp['dpT']*xsec*2*np.pi*2*0.8))
            jetscape[eloss][cent] = tmp
        ## construct the 0-20 and 20-40 centralities from the above
        tmp1 = util.combine_centralities(jetscape[eloss]['00-05'],jetscape[eloss]['05-10'])
        tmp2 = util.combine_centralities(tmp1, jetscape[eloss]['10-20'])
        jetscape[eloss]['00-20'] = tmp2
        tmp3 = util.combine_centralities(jetscape[eloss]['20-30'],jetscape[eloss]['30-40'])
        jetscape[eloss]['20-40'] = tmp3

    ## Read in the prompt and thermal calculations of JF Paquet
    channels = ['thermal','preEq','prompt']
    other_photons = {}
    tmpl_jf = "../../other_data/JF_MultiMessenger/PbPb2760_{cent}_{chan}.csv"
    for cent in ['00-20','20-40']:
        other_photons[cent] = {}
        for ch in channels:
            if ch == 'prompt':
                continue
            tmp = np.loadtxt(tmpl_jf.format(cent=cent,chan=ch),unpack=True,delimiter=',')
            x = tmp[0]
            y = tmp[1]
            other_photons[cent][ch] = pd.DataFrame({'pT':x,'N':y, 'dN':np.zeros_like(y)})
    kfactor = -1
    tmp = pd.read_csv('../../jetscape_data/prompt_photons/PbPb_2760/gamma_spectra.csv', comment='#')
    tmp['pT']   = 0.5*(tmp['pTmin'] + tmp['pTmax'])
    tmp = tmp[tmp['pT'].between(lower_bound['00-20'],upper_bound['00-20'])]
    tmp_dx  = tmp['pTmax'] - tmp['pTmin']
    tmp_y   = tmp['prompt']/(2*np.pi*tmp['pT']*tmp_dx*xsec*2*0.8)
    tmp_dy  = tmp['dprompt']/(2*np.pi*tmp['pT']*tmp_dx*xsec*2*0.8)
    for cent in ['00-20','20-40']:
        yy  = tmp_y*multiplicity[cent]
        dyy = tmp_dy*np.sqrt(multiplicity[cent])
        if kfactor < 0:
            f = interp1d(tmp['pT'], np.log(yy), kind='linear', fill_value='extrapolate')
            dat_x, dat_y = expdata['00-20']['x'].to_list(), expdata['00-20']['y'].to_list()
            kfactor = dat_y[-1]/np.exp(f(dat_x[-1]))
        print(kfactor)
        yy = kfactor*yy
        dyy = np.sqrt(kfactor)*dyy
        other_photons[cent]['prompt'] = pd.DataFrame({'pT':tmp['pT'].to_list(), 'N':yy.to_list(), 'dN':dyy.to_list()})
    
    x = expdata['00-20']['x']
    for cent in other_photons:
        specs = other_photons[cent]
        f_prompt  = interp1d(specs['prompt']['pT'] , np.log(specs['prompt']['N'])  , kind='linear', fill_value='extrapolate')
        f_thermal = interp1d(specs['thermal']['pT'], np.log(specs['thermal']['N']) , kind='linear', fill_value='extrapolate')
        f_preEq   = interp1d(specs['preEq']['pT']  , np.log(specs['preEq']['N'])   , kind='linear', fill_value='extrapolate')
        df_prompt = interp1d(specs['prompt']['pT'] , np.log(specs['prompt']['dN']) , kind='linear', fill_value='extrapolate')
        ytot      = np.exp(f_preEq(x)) + np.exp(f_prompt(x)) + np.exp(f_thermal(x)) 
        dytot     = np.exp(df_prompt(x))
        other_photons[cent]['total'] = (x,ytot,dytot)

    ratios = {}
    totals = {}
    for eloss in ['martini','cujet','none']:
        totals[eloss] = {}
        ratios[eloss] = {}
        for cent in ['00-20','20-40']:
            x,oth,doth         = other_photons[cent]['total']
            if eloss != 'none':
                tmpx  = jetscape[eloss][cent]['pT']
                spec  = jetscape[eloss][cent]['conv']+ jetscape[eloss][cent]['brem']
                dspec = np.sqrt(jetscape[eloss][cent]['dconv']**2 +jetscape[eloss][cent]['dbrem']**2)
                fjetmed     = interp1d(tmpx, np.log(spec), kind='linear', fill_value='extrapolate')
                dfjetmed    = interp1d(tmpx, np.log(dspec), kind='linear', fill_value='extrapolate')
                
                jet_medium_spec = np.exp(fjetmed(x))
                djet_medium_spec = np.exp(dfjetmed(x))
                total_spec  = oth + jet_medium_spec
                dtotal_spec = np.sqrt(doth**2 + djet_medium_spec**2)
                totals[eloss][cent] = (x, total_spec, dtotal_spec)

                ratio_j_t = jet_medium_spec/total_spec
                ratio_o_t = oth/total_spec
                ratios[eloss][cent] = (x, ratio_j_t, ratio_o_t)
            else:
                total_spec  = oth
                dtotal_spec = doth
                totals[eloss][cent] = (x, total_spec, dtotal_spec)
                ratio_j_t = 0/total_spec
                ratio_o_t = oth/total_spec
                ratios[eloss][cent] = (x, ratio_j_t, ratio_o_t)
    return ratios, totals

def take_ratio_to_data_and_plot(axis, data, jetscape, color, marker, zorder):
    spec, dspec = jetscape
    pT, dpT = data['x'], data['dx']
    ratio = spec / data['y']
    dratio_stat_plus  = ratio*np.sqrt(data['dy_stat+']**2/data['y']**2 + dspec**2/spec**2)
    dratio_stat_minus = ratio*np.sqrt(data['dy_stat+']**2/data['y']**2 + dspec**2/spec**2)
    dratio_syst_plus  = ratio*data['dy_syst+']/data['y']
    dratio_syst_minus = ratio*data['dy_syst-']/data['y']
    
    errorboxes_2 = [Rectangle((x-delx, y - abs(dy_minus)), width=2*delx, height=(abs(dy_minus)+abs(dy_plus)), zorder=zorder)
                for x, delx, y, dy_minus, dy_plus in zip(pT, 0.5*dpT, ratio, dratio_syst_minus, dratio_syst_plus)]
    pc2 = PatchCollection(errorboxes_2, facecolor=color, edgecolor=color, alpha=0.5)
    axis.add_collection(pc2)
    axis.errorbar(pT, ratio, xerr=0.5*dpT,yerr=[dratio_stat_minus, dratio_stat_plus], color=color,marker=marker,linewidth=1,fmt='none', zorder=zorder)
    axis.scatter(pT, ratio, color=color,marker=marker, zorder=zorder)