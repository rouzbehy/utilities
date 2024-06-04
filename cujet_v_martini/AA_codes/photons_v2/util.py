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
    "xtick.major.size" : 12,
    "ytick.major.size" : 12,
    "xtick.minor.size" : 6,
    "ytick.minor.size" : 6,
    "axes.spines.right": False,
    "axes.spines.top"  : False,
    "legend.frameon"   : False,
    "axes.labelsize" : 30}
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

## Plotting function for experimental data:
def plot_expr_data_on_axis(axis, data, marker,color=css['black'],face=css['dimgray'], s=40):
    deltax = 0.5*(data["xhigh"]-data["xlow"])
    errorboxes = [Rectangle((x-delx, y - yerrlow), width=2*delx, height=abs(yerrlow)+abs(yerrhigh), zorder=5)
                    for x, delx, y, yerrlow, yerrhigh in
                    zip(data["x"], deltax, data["y"], data["dy_syst+"], [-1*e if e < 0 else e for e in data["dy_syst-"]])]
     
    # Create patch collection with specified colour/alpha
    pc = PatchCollection(errorboxes, facecolor=face, alpha=0.2)
    axis.add_collection(pc)
    axis.scatter(data["x"], data["y"], color=color, marker=marker, s=s)
    axis.errorbar(data["x"], data["y"], xerr=deltax, 
                                        yerr=[[-1*e if e < 0 else e for e in data["dy_stat-"]], data["dy_stat+"]], 
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
    y = y[::2] + y[1::2]
    dy_raw = list(zip(dy[::2],dy[1::2]))
    dy = np.sqrt(np.array([sum(map(lambda x:x*x, e)) for e in dy_raw]))
    xmin = np.array(xmin)
    xmax = np.array(xmax)
    xmin_resummed = np.array([xmin[0], *xmax[1:-2:2]])
    xmax_resummed = np.array([*xmin[2::2], xmax[-1]])
    return xmin_resummed, xmax_resummed, y, dy

def combine_centralities(df1, df2):
    xmin = df1['ptmin']
    xmax = df1['ptmax']
    x = 0.5*(xmin+xmax)
    dx = xmax - xmin
    conv =0.5*(df1['conv']+ df2['conv']) 
    brem =0.5*(df1['brem']+ df2['brem']) 
    othr =0.5*(df1['othr']+ df2['othr']) 

    dconv =np.sqrt(df1['dconv']* df1['dconv']+ df2['dconv']* df2['dconv'])
    dbrem =np.sqrt(df1['dbrem']* df1['dbrem']+ df2['dbrem']* df2['dbrem'])
    dothr =np.sqrt(df1['dothr']* df1['dothr']+ df2['dothr']* df2['dothr'])
    return pd.DataFrame({'pT':x,'dpT':dx,'ptmin':xmin,'ptmax':xmax,
                          'conv':conv,'dconv':dconv,'brem':brem,'dbrem':dbrem,
                          'othr':othr,'dothr':dothr})
