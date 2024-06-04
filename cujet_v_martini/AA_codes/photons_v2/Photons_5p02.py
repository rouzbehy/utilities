import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.collections import PatchCollection
from matplotlib.patches import Rectangle
import util
import jetDicts as ddicts
from dictionaries import multiplicity_PbPb_5p02 as multiplicity
from dictionaries import oversampling_factors, channel_linestyles
from COLORS import module_colors, colors, channel_colors
import photonfuncs as pfuncs
my_rcParams={
    "text.usetex": True,
    "font.family": "Georgia",
    "font.size": 30,
    "lines.linewidth": 2,
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
plt.rcParams.update(my_rcParams)


if __name__=='__main__':

    sigma = 68 ## Inelastic cross section

    ## read in experimental data
    ROOT_DIRECTORY = "/Users/rmyazdi/Documents/research/jetscape_project/"

    exp_loc = ROOT_DIRECTORY + "expt/PbPb_5p02/Photons/fileDirPhotonSpectrum.txt"
    expdata = pd.read_csv(exp_loc,comment='#').rename(columns={
        'PT [GEV/c]':'x',
        'PTLOW':'xlow',
        'PTHIGH':'xhigh',
        '1/(2*PI*Nevt)*(d**2N)/(PT*dPT*dY) [(GEV/c)**-2]':'y',
        'stat':'dy_stat','sys':'dy_syst'})

    expdata['dx']       = expdata['xhigh']-expdata['xlow']
    expdata['dy_stat+'] = expdata['dy_stat']
    expdata['dy_stat-'] = expdata['dy_stat']
    expdata['dy_syst+'] = expdata['dy_syst']
    expdata['dy_syst-'] = expdata['dy_syst']
    expdata['dratio_stat_plus']  = expdata['dy_stat+']**2/expdata['y']**2 
    expdata['dratio_stat_minus'] = expdata['dy_stat+']**2/expdata['y']**2 
    expdata['dratio_syst_plus']  = expdata['dy_syst+']/expdata['y']
    expdata['dratio_syst_minus'] = expdata['dy_syst-']/expdata['y']

    tail1    = expdata.tail(1)
    x_vals_1 = expdata['x'] 
    last_data_point_1 = [tail1['x'].iloc[0], tail1['y'].iloc[0]]
    centralities = ['00-10','10-20','30-50']
    prompt, kompost, thermal, martini, cujet = {}, {}, {}, {}, {}
    kfactors = {}
    for cent in centralities:

        ## Read in the prompt photons
        prompt_loc = ROOT_DIRECTORY+"v2/jetscape_data/prompt_photons/PbPb_5020/"
        tmp_prompt, kfactor = pfuncs.get_prompt_photons(prompt_loc, multiplicity[cent], last_data_point_1, sigma, kfactor=1)
        #tmp_prompt = pfuncs.combine_centralities(tmp_prompt, tmp_prompt, 2, 'prompt')
        tmp_prompt = pfuncs.harmonize_x_axis(tmp_prompt, x_vals_1, 'prompt')
        prompt[cent] = tmp_prompt
        ## Thermal Photons:
        # THERMAL_ROOT_DIR = "/Users/rmyazdi/Documents/research/hydro_plots/VISHNU/HydroProfiles/"
        # tmp_thermal = pfuncs.get_thermal_photons('PbPb5020', cent, THERMAL_ROOT_DIR)
        # tmp_thermal = pfuncs.harmonize_x_axis(tmp_thermal, x_vals_1, channel='thermal')
        KOMPOST_ROOT_DIR = ROOT_DIRECTORY + 'v2/other_data/JF_MultiMessenger/'
        tmp_thermal = pfuncs.get_JF_photons('PbPb5020', cent, 'thermal', KOMPOST_ROOT_DIR)
        tmp_thermal = pfuncs.harmonize_x_axis(tmp_thermal, x_vals_1, channel='thermal')
        #tmp_thermal = pfuncs.combine_centralities(tmp_thermal, tmp_thermal, 2, 'thermal')
        #print(tmp_thermal.columns)
        thermal[cent] = tmp_thermal
        ## Pre-Equilibrium
        KOMPOST_ROOT_DIR = ROOT_DIRECTORY + 'v2/other_data/JF_MultiMessenger/'
        tmp_kompost = pfuncs.get_JF_photons('PbPb5020', cent, 'preEq', KOMPOST_ROOT_DIR)
        tmp_kompost = pfuncs.harmonize_x_axis(tmp_kompost, x_vals_1, channel='kompost')
        #tmp_kompost = pfuncs.combine_centralities(tmp_kompost, tmp_kompost, 2, 'kompost')
        #print(tmp_kompost.columns)
        kompost[cent] = tmp_kompost
        ## Jet Medium Photons:
        JMEDIUM_ROOT_DIR = ROOT_DIRECTORY + 'v2/jetscape_data/'
        tmp_martini   = pfuncs.get_jet_medium_photons('martini', cent, sigma, '5020', 'PbPb5020', multiplicity, JMEDIUM_ROOT_DIR)
        tmp_martini   = pfuncs.harmonize_x_axis(tmp_martini, x_vals_1, channel='jet-medium')
        #tmp_martini   = pfuncs.combine_centralities(tmp_martini, tmp_martini, 2, 'jet-medium')
        #print(tmp_martini.columns)
        martini[cent] = tmp_martini

        tmp_cujet = pfuncs.get_jet_medium_photons('cujet', cent, sigma, '5020', 'PbPb5020', multiplicity, JMEDIUM_ROOT_DIR)
        tmp_cujet = pfuncs.harmonize_x_axis(tmp_cujet, x_vals_1, channel='jet-medium')
        #tmp_cujet = pfuncs.combine_centralities(tmp_cujet, tmp_cujet, 2, 'jet-medium')
        #print(tmp_cujet.columns)
        cujet[cent] = tmp_cujet

    ## Construct the 0-20 % centrality spectra for each photon channel:
    cent1, cent2 = '00-10', '10-20'
    spectra = {'prompt':prompt, 'kompost':kompost, 'thermal':thermal, 'jet-medium':[martini, cujet]}
    for channel in spectra:
        if channel != 'jet-medium':
            spec = spectra[channel]
            spec['00-20'] = pfuncs.average_spectra(spec[cent1], spec[cent2], channel_names=[channel])
        else:
            spec_mart, spec_cujet = spectra[channel]
            spec_mart['00-20'] = pfuncs.average_spectra(spec_mart[cent1], spec_mart[cent2], channel_names=['conv','brem','total'])
            spec_cujet['00-20'] = pfuncs.average_spectra(spec_cujet[cent1], spec_cujet[cent2], channel_names=['conv','brem','total'])
        
    ## Construct totals:
    cent = '00-20'
    tail = prompt[cent].tail(1)
    k = last_data_point_1[1]/ tail['prompt'].iloc[0]
    print(f"kfactor: {k:0.5f}")
    prompt[cent]['prompt'] *= k
    totals = {}
    tot_no_jmed_y   = prompt[cent]['prompt'] + kompost[cent]['kompost'] + thermal[cent]['thermal'] 
    dtot_no_jmed_y  = prompt[cent]['dprompt']
    tot_martini_y   = tot_no_jmed_y + martini[cent]['total']
    dtot_martini_y  = np.sqrt(dtot_no_jmed_y** 2 + martini[cent]['dtotal']**2)
    tot_cujet_y     = tot_no_jmed_y + cujet[cent]['total']
    dtot_cujet_y    = np.sqrt(dtot_no_jmed_y** 2 + cujet[cent]['dtotal']**2)
    totals = pd.DataFrame({'pT':x_vals_1, 'No Jet-Med.':tot_no_jmed_y, 'err No Jet-Med.':dtot_no_jmed_y,
                                    'MATTER+MARTINI':tot_martini_y, 'err MATTER+MARTINI':dtot_martini_y,
                                    'MATTER+CUJET':tot_cujet_y, 'err MATTER+CUJET':dtot_cujet_y})

    ## Construct ratios:

    ratios = {}
    for cent in centralities:
        ratios[cent] = {}
        prompt_spectrum    = prompt[cent]['prompt']
        kompost_spectrum   = kompost[cent]['kompost']
        thermal_spectrum   = thermal[cent]['thermal']
        jet_medium_martini = martini[cent]['total']
        jet_medium_cujet   = cujet[cent]['total']
        total_with_martini = prompt_spectrum + thermal_spectrum + kompost_spectrum + jet_medium_martini
        total_with_cujet = prompt_spectrum + thermal_spectrum + kompost_spectrum + jet_medium_cujet

        prompt_ratio = prompt_spectrum/total_with_martini
        thermal_ratio = thermal_spectrum/total_with_martini
        kompost_ratio = kompost_spectrum/total_with_martini
        martini_ratio = jet_medium_martini/total_with_martini
        ratios[cent]['martini'] = pd.DataFrame({'pT':prompt[cent]['pT'], 'prompt':prompt_ratio,
                                                       'thermal':thermal_ratio,'kompost':kompost_ratio,
                                                       'jet-medium':martini_ratio}) 

        prompt_ratio = prompt_spectrum/total_with_cujet
        thermal_ratio = thermal_spectrum/total_with_cujet
        kompost_ratio = kompost_spectrum/total_with_cujet
        cujet_ratio = jet_medium_cujet/total_with_cujet
        ratios[cent]['cujet'] = pd.DataFrame({'pT':prompt[cent]['pT'], 'prompt':prompt_ratio,
                                                       'thermal':thermal_ratio,'kompost':kompost_ratio,
                                                       'jet-medium':cujet_ratio}) 
    ## Plotting time
    labels =['No Jet-Med.','MATTER+MARTINI','MATTER+CUJET']
    alpha=0.3
    fig1, (ax1, ax2)= plt.subplots(2, 1, height_ratios=(3,1), sharex=True, sharey='row')
    cent = r'$00$-$20$'
    data = expdata
    sims = totals
    ax1.set_yscale('log')
    ax1.text(0.05, 0.05, cent + r'$\%$', transform=ax1.transAxes)
    util.plot_expr_data_on_axis(ax1, data, marker='s', color='black', face='black', s=20)
    for label in labels: 
        col = module_colors[label]
        y, dy = sims[label], sims[f'err {label}']
        ax1.plot(sims['pT'], y, color=col, label=label)
        ax1.fill_between(sims['pT'], y+dy, y-dy, color=col, alpha=alpha)
        ratio = y/data['y']
        err = dy**2 / y**2
        dratio_plus_stat = ratio * np.sqrt(err+ data['dratio_stat_plus'])
        dratio_minus_stat = ratio * np.sqrt(err + data['dratio_stat_minus'])
        dratio_plus_syst = ratio * data['dratio_syst_plus']
        dratio_minus_syst = ratio * data['dratio_syst_minus']
        deltax = 0.5*(data["xhigh"]-data["xlow"])
        errorboxes = [Rectangle((x-delx, y - yerrlow), width=2*delx, 
                                height=abs(yerrlow)+abs(yerrhigh), zorder=5)
                        for x, delx, y, yerrlow, yerrhigh in
                        zip(data["x"], deltax, ratio, dratio_plus_syst, dratio_minus_syst)]
        pc = PatchCollection(errorboxes, facecolor=col, alpha=0.2)
        ax2.add_collection(pc)
        ax2.scatter(data["x"] , ratio, color=col, marker='p', s=10)
        ax2.errorbar(data["x"], ratio, xerr=deltax, yerr=[dratio_minus_stat, dratio_plus_stat], 
                                    fmt='none', 
                                    color=col)

    ax2.set_xlabel(r'$p_T$ (GeV)')
    ax2.axhline(y=1., linestyle='dotted', color='black')

    handles = [Line2D([],[],color=c, label=l) for l, c in module_colors.items() if l != 'MATTER']
    handles.append(Line2D([],[],linewidth=None, color='black', marker='s', markersize=5, label=r'ALICE (2022) $|\eta|<0.8$ [Preliminary]'))
    ax1.legend(loc='upper right', handles=handles, fontsize=25)
    ax2.set_ylabel('Theory over'+'\n data')
    ax1.set_ylabel(r'$\frac{1}{2\pi\;p_T}\frac{dN}{dp_T\;d\eta}$ (GeV)${}^{-2}$')

    ## Now plot channel ratios to total:
    fig2, axes2 = plt.subplots(3, 2, sharex=True, sharey=True)
    for irow, centrality in enumerate(ratios):
       specs = ratios[centrality]
       for icol, eloss in enumerate(specs):
           spec = specs[eloss]
           ax = axes2[irow][icol]
           if irow == 0:
               ax.set_title(eloss.upper())
           for channel in ['thermal','prompt','kompost']:
               ax.plot(spec['pT'], spec[channel], color=channel_colors[channel])
           ax.plot(spec['pT'], spec['jet-medium'], color=module_colors['MATTER+'+eloss.upper()])
           ax.text(0.05,0.8, centrality+r'$\%$', transform=ax.transAxes)
    for ax in axes2[2,:]:
       ax.set_xlabel(r'$p_T$ (GeV)')
    for i in range(3):
       axes2[i][0].set_ylabel(r'Channel to '+'\n'+ 'Total')
    handles = [Line2D([],[],color=channel_colors['prompt'], label='Prompt'),
              Line2D([],[],color=channel_colors['kompost'],label='Pre. Eq.'),
              Line2D([],[],color=channel_colors['thermal'],label='Thermal'),
              Line2D([],[],color=module_colors['MATTER+MARTINI'], label='MATTER+MARTINI'),
              Line2D([],[],color=module_colors['MATTER+CUJET'],label='MATTER+CUJET')]
    axes2[0][1].legend(loc='upper left', handles=handles, bbox_to_anchor=(1.05,0.9))
    plt.show()    
