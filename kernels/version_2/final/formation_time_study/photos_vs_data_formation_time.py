## import scipy and matplotlib modules:
import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.ticker as plticker
import matplotlib.pyplot as plt
from matplotlib.colors import CSS4_COLORS as css
from matplotlib.collections import PatchCollection
from matplotlib.lines import Line2D
from matplotlib.patches import Rectangle
#from scipy.interpolate import InterpolatedUnivariateSpline as interpolate
from scipy.interpolate import interp1d as interpolate
## my custom modules
import util
import dictionaries as my_dicts

plt.rcParams.update(util.my_rcParams)
from COLORS import rate_set_colors as rate_colors
from COLORS import colors
diff_xsec = 62.039
def get_spectrum(fname, oversampling_factor=1, centrality='00-05'):
    pTmin, pTmax = (0,19) if 'final' in fname else (0, 20)
    tmp = pd.read_csv(fname, comment='#')
    if 'pTmin' in tmp:
        tmp['pT'] = 0.5 * (tmp['pTmin'] + tmp['pTmax'])
        tmp['dpT'] = tmp['pTmax'] - tmp['pTmin']
    else:
        tmp['pT'] = 0.5 * (tmp['ptmin'] + tmp['ptmax'])
        tmp['dpT'] = tmp['ptmax'] - tmp['ptmin']
    tmp = tmp[tmp['pT'].between(pTmin, pTmax)]
    for col in ['conv','dconv','brem','dbrem']:
        tmp[col] /= (diff_xsec*oversampling_factor)
        tmp[col] *= my_dicts.multiplicity[centrality]
    tmp['total']  = tmp['conv'] + tmp['brem']
    tmp['dtotal'] = np.sqrt(tmp['dconv']**2 + tmp['dbrem']**2)
    tmp['total']  /= (2*np.pi*2*0.8*tmp['pT']*tmp['dpT'])
    tmp['dtotal'] /= (2*np.pi*2*0.8*tmp['pT']*tmp['dpT']) 
    tmp['conv']   /= (2*np.pi*2*0.8*tmp['pT']*tmp['dpT'])
    tmp['dconv']  /= (2*np.pi*2*0.8*tmp['pT']*tmp['dpT']) 
    tmp['brem']   /= (2*np.pi*2*0.8*tmp['pT']*tmp['dpT'])
    tmp['dbrem']  /= (2*np.pi*2*0.8*tmp['pT']*tmp['dpT'])  
    return tmp

## READ IN THE CALCULATION:
rates = [1,2,3]
cents = {'0_5':'00-05','5_10':'05-10','10_20':'10-20'}
template = '../martini_results/final_PbPb_2p76/rset_{r}/cent_{cent}/gamma_spectra.csv'
lo_nlo_np_work = {r: {cent: get_spectrum(template.format(r=r,cent=cent), 1000, cents[cent]) 
                        for cent in cents} for r in rates}

## Read in the jetscape work
prefix = '../../../../jetscape_project/v2/jetscape_data'
spec_dir = '/sqrt_s_2760/martini/PbPb_2760/'
fname_AA = 'PbPb2760_{cent}_photon_spec_0.80.csv'
template = prefix+spec_dir+fname_AA 
factor = {'00-05':1000000, '05-10':10000, '10-20':10000}
jetscape = {cent: get_spectrum(template.format(cent=cent), factor[cent]) for cent in factor} 

## Construct the 0-20 centrality line:
spectra_kernels = {}
for r in lo_nlo_np_work:
    tmp = lo_nlo_np_work[r]
    conv = (tmp['0_5']['conv'] + tmp['5_10']['conv'] + tmp['10_20']['conv'])/4
    dconv = np.sqrt(tmp['0_5']['dconv']*tmp['0_5']['dconv'] + tmp['5_10']['dconv']**2 + tmp['10_20']['dconv']**2)/4 
    brem = (tmp['0_5']['brem'] + tmp['5_10']['brem'] + tmp['10_20']['brem'])/4
    dbrem = np.sqrt(tmp['0_5']['dbrem']**2 + tmp['5_10']['dbrem']**2 + tmp['10_20']['dbrem']**2)/4 
    total = (tmp['0_5']['total'] + tmp['5_10']['total'] + tmp['10_20']['total'])/4
    dtotal = np.sqrt(tmp['0_5']['dtotal']**2 + tmp['5_10']['dtotal']**2 + tmp['10_20']['dtotal']**2)/4 
    x, dx = tmp['0_5']['pT'], tmp['0_5']['dpT']
    spectra_kernels[r] = pd.DataFrame({'pT':x.to_list(), 'dpT':dx.to_list(),
                                        'conv':conv.to_list(), 'dconv':dconv.to_list(),
                                        'brem':brem.to_list(), 'dbrem':dbrem.to_list(),
                                        'total':total.to_list(), 'dtotal':dtotal.to_list()})


total  = (jetscape['00-05']['total'] + jetscape['05-10']['total'] + jetscape['10-20']['total'])/4
dtotal = np.sqrt(jetscape['00-05']['dtotal']**2 + jetscape['05-10']['dtotal']**2 + jetscape['10-20']['dtotal']**2)/4

conv  = (jetscape['00-05']['conv'] + jetscape['05-10']['conv'] + jetscape['10-20']['conv'])/4
dconv = np.sqrt(jetscape['00-05']['dconv']**2 + jetscape['05-10']['dconv']**2 + jetscape['10-20']['dconv']**2)/4

brem  = (jetscape['00-05']['brem'] + jetscape['05-10']['brem'] + jetscape['10-20']['brem'])/4
dbrem = np.sqrt(jetscape['00-05']['dbrem']**2 + jetscape['05-10']['dbrem']**2 + jetscape['10-20']['dbrem']**2)/4

jetscape_martini = pd.DataFrame({'pT':jetscape['00-05']['pT'].to_list(), 
                                 'dpT':jetscape['00-05']['dpT'].to_list(),
                                 'conv':conv.to_list(),'dconv':dconv.to_list(),
                                 'brem':brem.to_list(),'dbrem':dbrem.to_list(),
                                 'total':total.to_list(),'dtotal':dtotal.to_list()})

## ALICE data:
exp_loc = "../../../exp_data/photons/HEPData-ins1394677-v1-Table_1.csv"
exp_data = pd.read_csv(exp_loc, comment='#')
exp_data = exp_data[exp_data['x'] > 2]
deltax = 0.5*(exp_data["xhigh"]-exp_data["xlow"])
#colours = my_dicts.rate_colours
ratenames = {1:'LO',2:'NLO',3:'NP'}

fig, axes = plt.subplots(2,1, height_ratios=[3,1], sharex=True, figsize=(16,9), gridspec_kw={'left':0.079, 'bottom':0.09,'right':0.995,'top':0.995, 'hspace':0.038})
util.plot_expr_data_on_axis(axes[0], exp_data, marker='s', color='black', face='black')

dexp_up = exp_data['dy_stat+']**2 + exp_data['dy_syst+']**2
dexp_down = exp_data['dy_stat-']**2 + exp_data['dy_syst-']**2
for rate in spectra_kernels:
    spec   = spectra_kernels[rate]
    colour = rate_colors[ratenames[rate]]
    fy  = interpolate(spec['pT'], np.log(spec['total']), kind='linear')
    dfy = interpolate(spec['pT'], np.log(spec['dtotal']), kind='linear')
    theo = np.exp(fy(exp_data['x']))
    dtheo = np.exp(dfy(exp_data['x']))

    ax = axes[0]
    ax.plot(exp_data['x'], theo, color=colour)
    ax.fill_between(exp_data['x'], theo + dtheo, theo-dtheo, color=colour, alpha=0.2)

    ax = axes[1]
    ratio = theo/exp_data['y']
    dratio_up = ratio*np.sqrt(dtheo**2/theo**2 + dexp_up/exp_data['y']**2)
    dratio_down = ratio*np.sqrt(dtheo**2/theo**2 + dexp_down/exp_data['y']**2)
     
    # Create patch collection with specified colour/alpha
    ax.scatter(exp_data["x"], ratio, color=colour, marker='*', s=60)
    ax.errorbar(exp_data["x"], ratio, xerr=deltax, 
                                        yerr=[dratio_down, dratio_up], 
                                        fmt='none', 
                                        color=colour)
    #ax.plot(exp_data['x'], ratio, color=colour)
    #ax.fill_between(exp_data['x'], ratio + dratio_up, ratio-dratio_down, color=colour, alpha=0.2)

# ax = axes[0]
# fy  = interpolate(jetscape_martini['pT'], np.log(jetscape_martini['total']), kind='linear')
# dfy = interpolate(jetscape_martini['pT'], np.log(jetscape_martini['dtotal']), kind='linear')
# theo = np.exp(fy(exp_data['x']))
# dtheo = np.exp(dfy(exp_data['x']))
# ax.plot(exp_data['x'], theo, color=rate_colors[6])
# ax.fill_between(exp_data['x'], theo + dtheo, theo-dtheo, color=colors[6], alpha=0.2)

# ax = axes[1]
# ratio = theo/exp_data['y']
# dratio_up = ratio*np.sqrt(dtheo**2/theo**2 + dexp_up/exp_data['y']**2)
# dratio_down = ratio*np.sqrt(dtheo**2/theo**2 + dexp_down/exp_data['y']**2)
# ax.scatter(exp_data["x"], ratio, color=colors[6], marker='*', s=60)
# ax.errorbar(exp_data["x"], ratio, xerr=deltax, yerr=[dratio_down, dratio_up], 
#                            fmt='none',color=colors[6])

axes[0].set_yscale('log')
axes[0].set_ylabel(r'$\frac{1}{2\pi p_T}\frac{dN^{\gamma}}{dp_T d\eta}$ GeV$^{-2}$')
axes[1].set_xlabel(r'$p^{\gamma}_T$ (GeV)')
axes[1].set_ylabel('Ratio to Data')


## Legend work:
handles = [Line2D([],[],color=c, label=f'Pythia+MARTINI({l})') for l, c in rate_colors.items()]
handles.append(Line2D([],[],color=colors[6], label='MATTER+MARTINI(LO)'))
#handles.append(Line2D([],[],label='ALICE (2016)', marker='s',color='black'))
axes[0].legend(loc='upper right', handles=handles)


## Read in the prompt, thermal and pre-equilibrium photons:
## Read in the prompt and thermal calculations of JF Paquet
channels = ['prompt','thermal','preEq']
nice_channels = {'prompt':'Prompt','thermal':'Thermal','preEq':'Pre-Equilibrium'}
jfp_spec = {}
tmpl_jf = "../../../../jetscape_project/v2/other_data/JF_MultiMessenger/PbPb2760_00-20_{chan}.csv"
for ch in channels:
    if ch == 'prompt':
        continue
    tmp = np.loadtxt(tmpl_jf.format(chan=ch),unpack=True,delimiter=',')
    x = tmp[0]
    y = tmp[1]
    jfp_spec[ch] = pd.DataFrame({'pT':x,'N':y})

tmp = pd.read_csv('../../../../jetscape_project/v2/jetscape_data/prompt_photons/PbPb_2760/gamma_spectra.csv', comment='#')
tmp_x   = 0.5*(tmp['pTmin'] + tmp['pTmax'])
tmp_dx  = tmp['pTmax'] - tmp['pTmin']
tmp_y   = tmp['prompt']/(2*np.pi*tmp_x*tmp_dx*diff_xsec*2*0.8)
tmp_dy  = tmp['dprompt']/(2*np.pi*tmp_x*tmp_dx*diff_xsec*2*0.8)
yy  = tmp_y*my_dicts.multiplicity['00-20']
dyy = tmp_dy*my_dicts.multiplicity['00-20']
## find k-factor:
#f = interpolate(tmp_x, np.log(yy), kind='linear', fill_value='extrapolate')
#dat_x, dat_y = exp_data['x'].to_list(), exp_data['y'].to_list()
#kfactor = dat_y[-1]/np.exp(f(dat_x[-1]))
#print(kfactor)
kfactor = 0.736
yy = yy * kfactor
dyy = dyy * np.sqrt(kfactor)
jfp_spec['prompt'] = pd.DataFrame({'pT':tmp_x, 'N':yy.to_list(), 'dN':dyy.to_list()})
## Rebin the JF results to match my bins (and the experimental data)
## as a sample. The x axis is shared among all of them:
rebinned_jf = {}
# pTmin, pTmax = exp_data['xlow'], exp_data['xhigh']
x = exp_data['x'].to_list()
for ch in channels:
    master_data = jfp_spec[ch]
    y = []
    f = interpolate(master_data['pT'],np.log(master_data['N']), kind='linear', fill_value='extrapolate') 
    for xval in x:
        yval = np.exp(f(xval))
        y.append(yval)
    
    new_master = pd.DataFrame({'pT':x,'N':y})
    rebinned_jf[ch] = new_master

total_non_jetmed = rebinned_jf['prompt']['N'] + rebinned_jf['thermal']['N'] + rebinned_jf['preEq']['N']

## piggybacking:
piggy_x = np.linspace(2,20,10)
piggy_jf = {}  
for ch in channels:
    master_data = jfp_spec[ch]
    y = []
    f = interpolate(master_data['pT'],np.log(master_data['N']), kind='linear', fill_value='extrapolate') 
    for xval in piggy_x:
        yval = np.exp(f(xval))
        y.append(yval)
    piggy_jf[ch] = np.array(y)
d_jf_reskinned = piggy_jf['prompt'] + piggy_jf['thermal'] + piggy_jf['preEq']

fig, axes = plt.subplots(2,1, height_ratios=[3,1], 
                        sharex=True, figsize=(16,9))
                        #gridspec_kw={'left':0.079, 'bottom':0.09,'right':0.995,'top':0.995, 'hspace':0.038})
util.plot_expr_data_on_axis(axes[0], exp_data, marker='s', color='black', face='black')

for rate in spectra_kernels:
    spec   = spectra_kernels[rate]
    colour = rate_colors[ratenames[rate]]
    fy  = interpolate(spec['pT'], np.log(spec['total']), kind='linear' , fill_value='extrapolate')
    dfy = interpolate(spec['pT'], np.log(spec['dtotal']), kind='linear', fill_value='extrapolate')
    
    theo  = np.exp(fy(exp_data['x']))
    theo = np.array([v1+v2 for v1,v2 in zip(theo, total_non_jetmed)])
    dtheo = np.exp(dfy(exp_data['x']))
    theo_1  = np.exp(fy(piggy_x))
    totaltt = theo_1 + d_jf_reskinned 
    # with open(f'./photon_sp/photons_00-20_rset_{rate}_2p76_martini_alone.csv', 'w') as f:
    #     f.write('pT,jmed,othr,total\n')
    #     for item in zip(piggy_x, theo_1, d_jf_reskinned, totaltt):
    #         line = [f'{v:0.5e}' for v in item]
    #         f.write(','.join(line)+'\n')
    ax = axes[0]
    ax.plot(exp_data['x'], theo, color=colour)
    ax.fill_between(exp_data['x'], theo + dtheo, theo-dtheo, color=colour, alpha=0.2)
    ax = axes[1]
    ratio       = theo/exp_data['y']
    dratio_up   = ratio*np.sqrt(dtheo**2/theo**2 + dexp_up/exp_data['y']**2)
    dratio_down = ratio*np.sqrt(dtheo**2/theo**2 + dexp_down/exp_data['y']**2)
    # Create patch collection with specified colour/alpha
    ax.scatter(exp_data["x"], ratio, color=colour, marker='*', s=60)
    ax.errorbar(exp_data["x"], ratio, xerr=deltax, 
                                        yerr=[dratio_down, dratio_up], 
                                        fmt='none', 
                                        color=colour)

# ax = axes[0]
# fy  = interpolate(jetscape_martini['pT'], np.log(jetscape_martini['total']), kind='linear')
# dfy = interpolate(jetscape_martini['pT'], np.log(jetscape_martini['dtotal']), kind='linear')
# theo = np.exp(fy(exp_data['x']))
# theo = np.array([v1+v2 for v1,v2 in zip(theo, total_non_jetmed)])
# dtheo = np.exp(dfy(exp_data['x']))
# ax.plot(exp_data['x'], theo, color=colors[6])
# ax.fill_between(exp_data['x'], theo + dtheo, theo-dtheo, color=colors[6], alpha=0.2)

# ax = axes[1]
# ratio = theo/exp_data['y']
# dratio_up = ratio*np.sqrt(dtheo**2/theo**2 + dexp_up/exp_data['y']**2)
# dratio_down = ratio*np.sqrt(dtheo**2/theo**2 + dexp_down/exp_data['y']**2)
# ax.scatter(exp_data["x"], ratio, color=colors[6], marker='*', s=60)
# ax.errorbar(exp_data["x"], ratio, xerr=deltax, yerr=[dratio_down, dratio_up], 
#                           fmt='none',color=colors[6])

axes[0].set_yscale('log')
axes[0].set_ylabel(r'$\frac{1}{2\pi p_T}\frac{dN^{\gamma}}{dp_T d\eta}$ GeV$^{-2}$')
axes[1].set_xlabel(r'$p^{\gamma}_T$ (GeV)')
axes[1].set_ylabel('Data over \n Theory' )


## Legend work:
handles = [Line2D([],[],color=c, label=f'Total (w. MARTINI({l}))') for l, c in rate_colors.items()]
#handles.append(Line2D([],[],color=rate_colors[6], label='Total (MATTER + MARTINI(LO))'))
handles.append(Line2D([],[],label='ALICE (2016)', marker='s',color='black'))
axes[0].legend(loc='upper right', handles=handles)
axes[0].text(0.05,0.30, r'Pb-Pb @ $\sqrt{s}=2.76$ ATeV', transform=axes[0].transAxes)
axes[0].text(0.05,0.2, r'0-20$\%$, $|\eta|<0.8$',transform=axes[0].transAxes)



# ## object names
# # total, dtotal, conv, dconv, brem, dbrem
# # total_non_jetmed, jfp_spec
# ## Ratios: jetmedium to total, channel-by-channel breakdown
# fig, axes = plt.subplots(1,2, sharex=True, figsize=(16,9), sharey=True,
#                           gridspec_kw={'left':0.079, 'bottom':0.09,\
#                                        'right':0.995,'top':0.995,'hspace':0.038})
# ax = axes[0]

# ## first: ratio of jet medium total channel to 

plt.show()
