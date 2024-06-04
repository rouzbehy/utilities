"""
    utility functions to declutter the main
    scripts.
    provide functionality for:
        - reading experimental files
        - reading in the theory calculation
        - performing integration if necessary (for v_2)
"""
from typing import Tuple
import pandas as pd
import numpy as np
from scipy.integrate import trapezoid, fixed_quad
from scipy.interpolate import InterpolatedUnivariateSpline as interpolator
from matplotlib.collections import PatchCollection
from matplotlib.patches import Rectangle
from matplotlib.colors import CSS4_COLORS as css
from dictionaries import exag_factor
## theory calculations:
my_rcParams={
    "text.usetex": True,
    "font.family": "Georgia",
    "font.size": 30,
    "lines.linewidth": 4,
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
def integrated_ratio(ax, data, theory_funcs, color, marker):
    """
        integrated_ratio(ax, data, theory):
            ax: axis to plot on
    """
    (ftheory, dftheory) = theory_funcs
    f = lambda x: np.exp(ftheory(np.log(x)))
    df = lambda x: np.exp(dftheory(np.log(x))) 
    area_theo = []
    area_theo_err = []
    data_dx =(data['xhigh'] - data['xlow'])
    data_dy_pos = np.sqrt(data['dy_stat+']*data['dy_stat+'] + data['dy_syst+']*data['dy_syst+'])
    data_dy_neg = np.sqrt(data['dy_stat-']*data['dy_stat-'] + data['dy_syst-']*data['dy_syst-']) 
    area_data = data_dx*data['y']
    area_data_high = data_dx*(data['y'] + data_dy_pos) - area_data
    area_data_low = area_data - data_dx*(data['y'] - data_dy_neg)
    ## now go bin by bin and calculate the area of the theory curve
    for i in range(len(data['y'])):
        area_theo.append(fixed_quad(f, a=data['xlow'].iat[i],b=data['xhigh'].iat[i], n=3)[0])
        area_theo_err.append(fixed_quad(df, a=data['xlow'].iat[i],b=data['xhigh'].iat[i], n=3)[0])
    
    area_theo = np.array(area_theo)
    area_theo_err = np.array(area_theo_err)
    ratio = area_theo/area_data
    #err_high = ratio*np.sqrt(area_data_low*area_data_low/(area_data*area_data)  + 0.25*area_theo_err*area_theo_err/(area_theo*area_theo))
    #err_low = ratio*np.sqrt(area_data_high*area_data_high/(area_data*area_data) + 0.25*area_theo_err*area_theo_err/(area_theo*area_theo)) 
    err_high = ratio*np.sqrt(0.25*area_theo_err*area_theo_err/(area_theo*area_theo))
    err_low =  ratio*np.sqrt(0.25*area_theo_err*area_theo_err/(area_theo*area_theo)) 


    #print("pT \t Data high \t Data low \t Theory \t Propagated High \t Propagated Low")
    #for item in zip(data['x'],area_data_high/area_data, area_data_low/area_data, area_theo_err, err_high, err_low):
    #    print("{:0.2f} \t {:0.3e} \t {:0.3e} \t {:0.3e} \t {:0.3e} \t {:0.3e}".format(*item))


    #err_low = ratio - area_theo/area_data_high
    #err_high =   area_theo/area_data_low - ratio
    ax.scatter(data['x'], ratio, color=color, marker=marker, s=45)
    ax.errorbar(data['x'], ratio, xerr=0.5*data_dx, 
                                 yerr=[err_low, err_high], 
                                 fmt='none', 
                                 color=color)

def integrated_ratio_v2(ax, data, theory_funcs, color):
    """
        integrated_ratio(ax, data, theory):
            ax: axis to plot on
    """
    ftheory = theory_funcs
    f = lambda x: np.exp(ftheory(np.log(x)))
    area_theo = []
    data_dx =(data['xhigh'] - data['xlow'])
    area_data = data_dx*data['y']
    for i in range(len(data['y'])):
        area_theo.append(fixed_quad(f, a=data['xlow'].iat[i],b=data['xhigh'].iat[i], n=3)[0])
    ratio = area_theo/area_data
    ax.plot(data['x'], ratio, color=color)

def get_photon_spectra_jetscape(eloss : str, centrality : str) -> dict:
    """
        get_photon_spectra_jetscape(eloss, centrality):
            eloss: energy loss module used to gen the data
            centrality: centrality class of the desired data
            returns a datafram
    """
    xSec  = 62.039
    fname = f"../data/{eloss}/{eloss}_{centrality}/photon_spec.csv"
    dat = pd.read_csv(fname, comment='#')
    oversampling_factor = exag_factor[eloss]
    pT    = 0.5*(dat['pTmin'] + dat['pTmax'])
    dpT   = dat['pTmax']-dat['pTmin']
    conv  = dat['conv'] /(2*np.pi*(2*1)*pT*dpT*oversampling_factor*xSec)
    dconv = dat['dconv']/(2*np.pi*(2*1)*pT*dpT*oversampling_factor*xSec)
    brem  = dat['brem'] /(2*np.pi*(2*1)*pT*dpT*oversampling_factor*xSec)
    dbrem = dat['dbrem']/(2*np.pi*(2*1)*pT*dpT*oversampling_factor*xSec)
    return {'pT': pT ,'conv':conv,'dconv':dconv,'brem':brem,'dbrem':dbrem}

def get_photon_v2_jetscape(eloss : str, centrality : str) -> dict:
    """
        get_photon_v2_jetscape(eloss, centrality):
            eloss: energy loss module used to gen the data
            centrality: centrality class of the desired data
            returns a tuple of spectra
    """
    fname = f"../data/{eloss}_{centrality}/photon_v2.csv"
    dat = pd.read_csv(fname, comment='#')
    pT = 0.5*(dat['pTmin'] + dat['pTmax'])
    dpT = dat['pTmax']-dat['pTmin']
    f_conv_c2 = interpolator(dat['pT'], dat['conv_c2']/(2*np.pi*(2*1)*pT*dpT),k=1,ext=0) ## 1: eta window
    f_conv_s2 = interpolator(dat['pT'], dat['conv_s2']/(2*np.pi*(2*1)*pT*dpT),k=1,ext=0) ## 1: eta window
    f_brem_c2 = interpolator(dat['pT'], dat['brem_c2']/(2*np.pi*(2*1)*pT*dpT),k=1,ext=0) ## 1: eta window
    f_brem_s2 = interpolator(dat['pT'], dat['brem_s2']/(2*np.pi*(2*1)*pT*dpT),k=1,ext=0) ## 1: eta window

    df_conv_c2 = interpolator(dat['pT'], dat['dconv_c2']/(2*np.pi*(2*1)*pT*dpT),k=1,ext=0) ## 1: eta window
    df_conv_s2 = interpolator(dat['pT'], dat['dconv_s2']/(2*np.pi*(2*1)*pT*dpT),k=1,ext=0) ## 1: eta window
    df_brem_c2 = interpolator(dat['pT'], dat['dbrem_c2']/(2*np.pi*(2*1)*pT*dpT),k=1,ext=0) ## 1: eta window
    df_brem_s2 = interpolator(dat['pT'], dat['dbrem_s2']/(2*np.pi*(2*1)*pT*dpT),k=1,ext=0) ## 1: eta window

    return {'pT': dat['pT'],'conv_c2':f_conv_c2,'dconv_c2':df_conv_c2,
            'conv_s2':f_conv_s2,'dconv_s2':df_conv_s2,
            'brem_c2':f_brem_c2,'dbrem_c2':df_brem_c2,
            'brem_s2':f_brem_s2,'dbrem_s2':df_brem_s2}
                                         
    return dat

def get_jet_spectra(eloss : str, centrality : str) -> Tuple:
    """
        get_jet_spectra(eloss, centrality):
            eloss: energy loss module used to gen the data
            centrality: centrality class of the desired data
            returns a tuple of spectrum functions
    """
    fname = f"../data/{eloss}/{eloss}_{centrality}/jets.csv"
    dat = pd.read_csv(fname, comment='#')
    return dat

def get_charged_spectra(eloss : str, centrality : str) -> dict:
    """
        get_charged_spectra(eloss, centrality):
            eloss: energy loss module used to gen the data
            centrality: centrality class of the desired data
            returns a dictionary of the interpolating 
            functiongs N and dN (N = differential yield)
    """
    xSec  = 62.039
    fname = f"../data/{eloss}/{eloss}_{centrality}/charged.csv"
    dat   = pd.read_csv(fname, comment='#')
    #f_ch  = interpolator(dat['pT'] , dat['N']/(2*np.pi*(2*1)*dat['pT']*dat['dpT']*xSec), k=1, ext=0)
    #df_ch = interpolator(dat['pT'], dat['dN']/(2*np.pi*(2*1)*dat['pT']*dat['dpT']*xSec),k=1,ext=0) 
    ##TODO: change the eta to 0.8 for the new data set
    dpT = dat['pTmax'] - dat['pTmin']
    pT = 0.5*(dat['pTmax'] + dat['pTmin'])
    Nch  = dat['N'] /(2*np.pi*(2*0.8)*pT*dpT*xSec)
    dNch = dat['dN']/(2*np.pi*(2*0.8)*pT*dpT*xSec)
    #return {'pT': dat['pT'],'N':f_ch, 'dN':df_ch}
    return {'pT':pT, 'N':Nch, 'dN':dNch}

def get_num_photons_jetscape(eloss: str, centrality : str) -> dict:
    """
        get_num_photons_jetscape(eloss, centrality):
            eloss: energy loss module used to gen the data
            centrality: centrality class of the desired data
            returns a 
    """
    photons = get_photon_spectra_jetscape(eloss, centrality)
    (n_conv, dn_conv, n_brem, dn_brem) = (1, 0, 1, 0)
    if eloss != 'cujet':
        n_conv = trapezoid(photons['conv'], x=photons['pT'])
        dn_conv =  trapezoid(photons['dconv'], x=photons['pT'])
        n_brem = trapezoid(photons['brem'], x=photons['pT'])
        dn_brem = trapezoid(photons['brem'], x=photons['pT'])
    else:
        n_conv = trapezoid(photons['conv'], x=photons['pT'])
        dn_conv =  trapezoid(photons['dconv'], x=photons['pT'])

    return {"conv":(n_conv,dn_conv),"brem":(n_brem,dn_brem)}


def combine_systematics(data):
    tmp1 = [np.sqrt(x*x + y*y) for (x,y) in zip(data['dy_syst1+'],data['dy_syst2+'])]
    tmp2 = [np.sqrt(x*x + y*y) for (x,y) in zip(data['dy_syst1-'],data['dy_syst2-'])]
    data['dy_syst+'] = tmp1
    data['dy_syst-'] = tmp2
    return data

## Plotting function for experimental data:
def plot_expr_data_on_axis(axis, data, marker,color=css['black'],face=css['dimgray'], s=40, factor=1):
    deltax = 0.5*(data["xhigh"]-data["xlow"])
    errorboxes = [Rectangle((x-delx, y - yerrlow), width=2*delx, height=abs(yerrlow)+abs(yerrhigh))
                    for x, delx, y, yerrlow, yerrhigh in
                    zip(data["x"], deltax, factor*data["y"], factor*abs(data["dy_syst+"]), factor*abs(data["dy_syst-"]))]
     
    # Create patch collection with specified colour/alpha
    pc = PatchCollection(errorboxes, facecolor=face, alpha=0.2)
    axis.add_collection(pc)
    axis.scatter(data["x"], factor*data["y"], color=color, marker=marker, s=s)
    axis.errorbar(data["x"], factor*data["y"], xerr=deltax, 
                                        yerr=[abs(data["dy_stat-"]), data["dy_stat+"]], 
                                        fmt='none', 
                                        color=color)

def plot_expr_v2_on_axis(axis, data, marker):
    deltax = 0.5*(data["xhigh"]-data["xlow"])
    axis.scatter(data["x"], data["y"], color=css['black'], marker=marker, s=55)
    axis.errorbar(data["x"], data["y"], xerr=deltax, 
                                        yerr=[-1*data["dy_tot-"], data["dy_tot+"]], 
                                        fmt='none', 
                                        color=css['black'])
def combine_bins(xmin, xmax, y, dy):
    y = np.array(y)
    nbins = 4
    y = y[::nbins] + y[nbins-1::nbins]
    dy_raw = list(zip(dy[::nbins],dy[nbins-1::nbins]))
    dy = np.sqrt(np.array([sum(map(lambda x:x*x, e)) for e in dy_raw]))
    xmin = np.array(xmin)
    xmax = np.array(xmax)
    xmin_resummed = np.array([xmin[0], *xmax[nbins-1:-nbins:nbins]])
    xmax_resummed = np.array([*xmin[nbins::nbins], xmax[-nbins]])
    print(xmin)
    print(xmin_resummed)
    print(xmax)
    print(xmax_resummed)
    return xmin_resummed, xmax_resummed, y, dy
