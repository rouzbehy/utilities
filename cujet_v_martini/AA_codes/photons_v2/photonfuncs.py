#!/usr/bin/env python3
import pandas as pd
import numpy as np
from scipy.interpolate import InterpolatedUnivariateSpline
from dictionaries import oversampling_factors

def get_jet_medium_photons(eloss, cent, sigma, energy, collision_system, multiplicity, ROOT_DIRECTORY='./'):
    tag = ''
    eta = 0
    if collision_system == 'PbPb2760':
        tag = 'PbPb_2760'
        eta = 0.8
    elif collision_system == 'PbPb5020':
        tag = 'PbPb_5020'
        eta = 0.8
    else:
        tag = 'AuAu_200'
        eta = 0.35
    fname_tmpl = ''
    if collision_system == 'PbPb2760':
        fname_tmpl = ROOT_DIRECTORY + 'sqrt_s_{energy}/{eloss}/{tag}/{coll}_{cent}_photon_spec_0.80.csv'
    elif collision_system == 'AuAu200':
        fname_tmpl = ROOT_DIRECTORY + 'sqrt_s_{energy}/AuAu/{eloss}/{cent}/photon_spec_0.35.csv'
    else:
        fname_tmpl = ROOT_DIRECTORY + 'sqrt_s_{energy}/maxt_200/{eloss}/{coll}_{cent}_photon_spec_0.80.csv' 

    fname = fname_tmpl.format(eloss=eloss,cent=cent, energy=energy, tag=tag, coll=collision_system)
    factor = 1
    if '2760' in collision_system:
        factor = oversampling_factors[f'{eloss}-{cent}']
    else:
        factor = 2000

    Nbin = multiplicity[cent]
    tmp = pd.read_csv(fname,comment='#')
    xmin = tmp['ptmin']
    xmax = tmp['ptmax']
    x = 0.5*(xmin+xmax)
    dx = xmax - xmin
    tmp['conv']  *= Nbin/(factor*(x*dx*sigma*2*np.pi*2*eta))
    tmp['dconv'] *= np.sqrt(Nbin)/(factor*(x*dx*sigma*2*np.pi*2*eta))
    tmp['brem']  *= Nbin/(factor*(x*dx*sigma*2*np.pi*2*eta))
    tmp['dbrem'] *= np.sqrt(Nbin)/(factor*(x*dx*sigma*2*np.pi*2*eta))
    tmp['total'] = tmp['conv'] + tmp['brem']
    tmp['dtotal'] = np.sqrt(tmp['dconv']**2 + tmp['dbrem']**2)
    tmp['pT'] = x
    tmp['dpT'] = dx
    return tmp

def combine_centralities(cent1, cent2, factor, channel):
    new_dataframe = {'pT':cent1['pT']}
    if channel == 'jet-medium':
        for ch in ['conv', 'brem', 'total']:
            y = (cent1[ch] + cent2[ch])/factor
            dy = np.sqrt(cent1[f'd{ch}']**2 + cent2[f'd{ch}']**2)/factor
            new_dataframe[ch] = y
            new_dataframe[f'd{ch}'] = dy
    elif channel == 'prompt':
        ch = 'N'
        y = (cent1[ch] + cent2[ch])/factor
        dy = np.sqrt(cent1[f'd{ch}']**2 + cent2[f'd{ch}']**2)/factor
        new_dataframe[ch] = y
        new_dataframe[f'd{ch}'] = dy
    else:
        ch = 'N'
        y = (cent1[ch] + cent2[ch])/factor
        new_dataframe[ch] = y
    return pd.DataFrame(new_dataframe)

def get_thermal_photons(coll_syst, cent, ROOT_DIRECTORY='./'):
    tmpl = ROOT_DIRECTORY + '{coll_syst}-Avg-{cent}/event-0/photon_total_Spvn.dat'
    tmp   = pd.read_csv(tmpl.format(coll_syst=coll_syst, cent=cent), dtype=np.float64)
    
    tmp.columns   = ['pT','N','tmp']
    dx1 = np.array([0.2])
    dx = np.diff(tmp['pT'])
    dx1= np.concatenate((dx1, dx))
    tmp['dpT'] = dx1
    tmp.fillna(0, inplace=True)
    tmp.drop('tmp', inplace=True, axis=1)
    return tmp

def get_JF_photons(coll_syst, cent, channel, ROOT_DIRECTORY='./'):
    tmpl_name = ROOT_DIRECTORY + "{coll_syst}_{cent}_{chan}.csv"
    tmp = np.loadtxt(tmpl_name.format(cent=cent,chan=channel, coll_syst=coll_syst),unpack=True,delimiter=',')
    x = tmp[0]
    y = tmp[1]
    return pd.DataFrame({'pT':x, 'N':y})

def harmonize_x_axis(dframe, new_xvals, channel):
    result = {}
    x = dframe['pT']
    result['pT'] = new_xvals
    if channel == 'jet-medium':
        tmp = dframe[dframe['total'] > 0]
        x = tmp['pT']
        for ch in ['conv', 'brem', 'total']:
            f = InterpolatedUnivariateSpline(x, np.log(tmp[ch]), k=4, ext='extrapolate')
            df = InterpolatedUnivariateSpline(x, np.log(tmp[f'd{ch}']), k=4, ext='extrapolate')
            y  = np.exp(f(new_xvals))
            dy = np.exp(df(new_xvals))
            result[ch] = y
            result[f'd{ch}'] = dy
    else:
         f = InterpolatedUnivariateSpline(x, np.log(dframe['N']), k=3, ext='extrapolate')
         y = np.exp(f(new_xvals))
         result[channel] = y
         if channel == 'prompt':
            df = InterpolatedUnivariateSpline(x, np.log(dframe[f'dN']), k=3, ext='extrapolate')
            dy = np.exp(df(new_xvals))
            result[f'd{channel}'] = dy
    return pd.DataFrame(result)

def get_prompt_photons(directory, Nbin, last_data_point, xsec, kfactor=None):
    tmp     = pd.read_csv(directory+'gamma_spectra.csv', comment='#')
    tmp     = tmp[tmp['prompt'] > 0]
    eta     = 0.8 if 'PbPb' in directory else 0.35
    tmp_x   = 0.5*(tmp['pTmin'] + tmp['pTmax'])
    tmp_dx  = tmp['pTmax'] - tmp['pTmin']
    tmp_y   = tmp['prompt']/(2*np.pi*tmp_x*tmp_dx*xsec *2*eta)
    tmp_dy  = tmp['dprompt']/(2*np.pi*tmp_x*tmp_dx*xsec*2*eta)
    yy  = tmp_y*Nbin
    dyy = tmp_dy*np.sqrt(Nbin)
    # find k-factor:
    k = kfactor
    if not kfactor:
        f = InterpolatedUnivariateSpline(tmp_x, np.log(yy), k=3, ext='extrapolate')
        dat_x, dat_y = last_data_point
        k = dat_y/np.exp(f(dat_x))
    yy = yy * k
    dyy = dyy * np.sqrt(k)
    return pd.DataFrame({'pT':tmp_x, 'dpT':tmp_dx, 'N':yy.to_list(), 'dN':dyy.to_list()}), k


def average_RHIC_data(dframe1, dframe2):
    """
        I've already massaged the PHENIX data by equalizing the x axes. so this is not 
        a very general function. I'm also hard coding the column names.
    """
    x = dframe1['$p_T$']
    xlow = dframe1['$p_T$ LOW']
    xhigh = dframe1['$p_T$ HIGH']
    
    y = 0.5 * (dframe1['inv.yield'] + dframe2['inv.yield'])
    tag = '$p_T$-uncorrelated-error +'
    uncorr_err_pos = 0.5*np.sqrt(dframe1[tag]**2 + dframe2[tag]**2)
    tag = '$p_T$-uncorrelated-error -'
    uncorr_err_neg = 0.5*np.sqrt(dframe1[tag]**2 + dframe2[tag]**2)
    tag = '$p_T$-correlated-error +'
    cor_err_pos = 0.5*np.sqrt(dframe1[tag]**2 + dframe2[tag]**2)
    tag = '$p_T$-correlated-error -'
    corr_err_neg = 0.5*np.sqrt(dframe1[tag]**2 + dframe2[tag]**2)

    return pd.DataFrame({'x':x, 'xlow':xlow, 'xhigh':xhigh, 'y':y, 
                         'dy_syst+':cor_err_pos, 'dy_syst-':corr_err_neg, 
                         'dy_stat+':uncorr_err_pos, 'dy_stat-':uncorr_err_neg})

def average_spectra(dframe1, dframe2, channel_names=''):
    result = {}
    result['pT'] = dframe1['pT']
    for channel_name in channel_names:
        result[channel_name] = 0.5*(dframe1[channel_name] + dframe2[channel_name])
        dname = f'd{channel_name}'
        if dname in dframe1:
            result[dname] = np.sqrt(dframe1[dname]**2 + dframe2[dname]**2)
    return pd.DataFrame(result)