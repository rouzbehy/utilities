import pandas as pd
from numpy import sqrt, pi, log
from scipy.interpolate import interp1d
from matplotlib.collections import PatchCollection
from matplotlib.patches import Rectangle

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

xsecs = {200: 42., 2760: 62, 5020:68}
kfactors = {200:1.12, 2760:0.736, 5020:0.95}

jet_medium_locations = {200: '../../jetscape_data/sqrt_s_200/{nuclei}/{eloss}/{cent}/photon_spec_0.35.csv',
                       2760: '../../jetscape_data/sqrt_s_2760/{eloss}/{syst_with_uscore}/{syst}_{cent}_photon_spec_0.80.csv',
                       5020: '../../jetscape_data/sqrt_s_5020/maxt_200/{eloss}/{syst}_{cent}_photon_spec_0.80.csv'}

prompt_locations = {200: '../../jetscape_data/prompt_photons/AuAu_200/gamma_spectra.csv',
                   2760: '../../jetscape_data/prompt_photons/PbPb_2760/gamma_spectra.csv',
                   5020: '../../jetscape_data/prompt_photons/PbPb_5020/gamma_spectra.csv'}

binaryColl_PbPb_2p76 = {'00-05':1615,
                '05-10':1268  , '10-15':1002,
                '15-20':790.1 , '20-30':548.5,
                '30-40':320.9 , '40-50':174.85,
                '10-20':896.05, '00-20':1168.77}

binaryColl_PbPb_5p02 = {'00-05':1762,'05-10':1380,
                          '10-15':1088,'00-10':(1762+1380)/2.,
                          '15-20':855.3,'00-20':(1762+1380+1088+855.3)/4.,
                          '10-20':(1088+855.3)/2,'30-50':(392.9+294.5+216.4+155.5)/4}

binaryColl_AuAu_200 = {'00-10':942.2, '10-20':591.55,'00-05':1053,'05-10':831.4}
binary = {200:binaryColl_AuAu_200, 2760:binaryColl_PbPb_2p76, 5020:binaryColl_PbPb_5p02}

## Oversampling factors:
oversample_PbPb_2p76 = {'martini-00-05':1000000,'martini-20-30':1000000,
            'martini-30-40':1000000 ,'martini-05-10':10000,
            'martini-10-20':10000   ,'martini-40-50':10000,
            'cujet-00-05':5000,'cujet-05-10':5000,
            'cujet-10-20':5000,'cujet-20-30':5000,
            'cujet-30-40':5000,'cujet-40-50':5000}



path ='../../other_data/JF_Photon_Calc/Spectra/'
path += 'thermal_from_{source}/{system}/C{cent}/average_sp.dat'

pTmin, pTmax = 1, 30
def construct_photon_spec(nuclei, energy, centrality, eloss, ptcut):
    """
        def construct_photon_spec(nuclei, energy, centrality, eloss):
            format of the input:
            *nuclei: string, AA
            *energy: integer
            *centrality: xx-yy
            *eloss: martini, cujet, none
    """
    ## Channels that are always there:
    eta = 0.35 if energy ==200 else 0.8
    prompt = pd.read_csv(prompt_locations[energy],comment='#')
    prompt['pT'] = 0.5* (prompt['pTmin']+prompt['pTmax'])
    factor = binary[energy][centrality]*kfactors[energy]
    prompt['prompt'] *= factor/xsecs[energy]
    prompt['dprompt'] *= sqrt(factor)/xsecs[energy]
    prompt = prompt[prompt['pT'].between(pTmin,pTmax)]
    prompt['dpT'] = prompt['pTmax'] - prompt['pTmin']
    prompt['N'] = prompt['prompt']/(2*pi*prompt['dpT']*prompt['pT']*2*eta)
    prompt['dN'] = prompt['dprompt']/(2*pi*prompt['dpT']*prompt['pT']*2*eta)

    tmp = [int(x) for x in centrality.split('-')]
    cent_for_jf = f'{tmp[0]}-{tmp[1]}'
    thermal_file = path.format(source='hydro', system=nuclei+f'{energy}', cent=cent_for_jf)
    thermal = pd.read_csv(thermal_file, header=None, usecols=[0,1], delimiter=r"\s+", names=['pT','N'])
    thermal.fillna(0.0, inplace=True)
    preeq_file = path.format(source='pre_equilibrium', system=nuclei+f'{energy}', cent=cent_for_jf)
    pre_Eq = pd.read_csv(preeq_file, header=None, usecols=[0,1], delimiter=r"\s+", names=['pT','N'])
    pre_Eq.fillna(0.0, inplace=True)
    ## determine what channels are to be included:
    jet_medium = 0
    if eloss != 'none':
        fname = jet_medium_locations[energy]
        if energy == 200:
            fname = fname.format(nuclei=nuclei, eloss=eloss, cent=centrality)
        elif energy == 2760:
            syst_with_uscore = nuclei+'_'+f'{energy}'
            fname = fname.format(eloss=eloss, syst_with_uscore=syst_with_uscore,
                                syst=nuclei+f'{energy}', cent=centrality)
        else:
            fname = fname.format(eloss=eloss, syst=nuclei+f'{energy}', cent=centrality) 
        tmp = pd.read_csv(fname, comment='#')
        tmp['pT'] = 0.5*(tmp['ptmax']+tmp['ptmin'])
        tmp = tmp[tmp['pT'].between(pTmin, pTmax)]
        tmp['dpT'] = tmp['ptmax'] - tmp['ptmin']
        conv, dconv = tmp['conv'], tmp['dconv']
        brem, dbrem = tmp['brem'], tmp['dbrem']
        oversampling = 2000 if energy != 2760 else oversample_PbPb_2p76[f'{eloss}-{centrality}']
        factor = binary[energy][centrality]/(xsecs[energy]*oversampling)
        conv = conv*factor/(2*pi*2*eta*tmp['pT']*tmp['dpT'])
        dconv = sqrt(factor)*dconv/(2*pi*2*eta*tmp['pT']*tmp['dpT'])
        brem = brem * factor /(2*pi*2*eta*tmp['pT']*tmp['dpT'])
        dbrem = dbrem * sqrt(factor)/(2*pi*2*eta*tmp['pT']*tmp['dpT'])
        tot = conv + brem
        dtot = sqrt(dconv**2 + dbrem**2)
        jet_medium = pd.DataFrame({'pT':tmp['pT'].to_list(),
                                 'dpT' :tmp['dpT'].to_list(),
                                   'N' :tot, 'dN':dtot})

    ## Now integrate the various spectra, no need to interpolate,
    ## just place cuts pTmin = [1, 5] for the two plots of PHENIX
    value1, value5 = 0, 0
    thermal1 = thermal[thermal['pT'] > ptcut]
    thermal_dNdy = 2*pi*sum(thermal1['pT']*0.2*thermal1['N'])
    preEq1 = pre_Eq[pre_Eq['pT'] > ptcut]
    preEq_dNdy = 2*pi*sum(preEq1['pT']*0.2*preEq1['N'])
    prompt1 = prompt[prompt['pT'] > ptcut]
    prompt_dNdy = 2*pi*sum(prompt1['pT']*prompt1['dpT']*prompt1['N'])
    dprompt_dNdy = 2*pi*sum(prompt1['pT']*prompt1['dpT']*prompt1['dN'])
    values1 = thermal_dNdy + preEq_dNdy + prompt_dNdy
    dvalues1 = dprompt_dNdy
    if eloss != 'none':
        jet_medium1 = jet_medium[jet_medium['pT'] > ptcut]
        values1 += 2*pi*sum(jet_medium1['pT']*jet_medium1['dpT']*jet_medium1['N'])
        dtmp = 2*pi*sum(jet_medium1['pT']*jet_medium1['dpT']*jet_medium1['dN'])
        dvalues1 = sqrt(dvalues1**2 + dtmp**2)
     
    return values1, dvalues1


def plot_with_error_box(axis, xvals, yvals, dxvals, dyvals, color):
    errorboxes = [Rectangle((x-delx, y - yerrlow), width=2*delx, height=abs(yerrlow)+abs(yerrhigh), zorder=5)
                    for x, delx, y, yerrlow, yerrhigh in
                    zip(xvals, dxvals, yvals, dyvals, -1*[v for v in dyvals])]
     
    # Create patch collection with specified colour/alpha
    pc = PatchCollection(errorboxes, facecolor=color, alpha=0.5)
    axis.add_collection(pc)
