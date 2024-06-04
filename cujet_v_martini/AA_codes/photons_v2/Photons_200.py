import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.collections import PatchCollection
from matplotlib.patches import Rectangle
import util
import jetDicts as ddicts
from dictionaries import multiplicity_AuAu_200 as multiplicity
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

    sigma = 42 ## Inelastic cross section

    ## read in experimental data
    ROOT_DIRECTORY = "/Users/rmyazdi/Documents/research/jetscape_project/"
    exp_loc = ROOT_DIRECTORY + "expt/AuAu_200/PHENIX_Photons/HEPData-ins1116179-v1-csv/PHENIX_{cent}.csv"
    centralities = ['00-05', '05-10', '10-15', '15-20']
    expt_data = {cent: pd.read_csv(exp_loc.format(cent=cent), comment='#') for cent in centralities}

    dat_00_10 = pfuncs.average_RHIC_data(expt_data['00-05'], expt_data['05-10'])
    dat_10_20 = pfuncs.average_RHIC_data(expt_data['10-15'], expt_data['15-20'])
    
    expdata = {'00-10': dat_00_10, '10-20':dat_10_20}

    for cent in expdata:
        tmp = expdata[cent]
        tmp['dx'] = tmp['xhigh']-tmp['xlow']
        tmp['dratio_stat_plus'] = tmp['dy_stat+']**2/tmp['y']**2 
        tmp['dratio_stat_minus'] = tmp['dy_stat+']**2/tmp['y']**2 
        tmp['dratio_syst_plus'] = tmp['dy_syst+']/tmp['y']
        tmp['dratio_syst_minus'] = tmp['dy_syst-']/tmp['y']
        expdata[cent] = tmp

    tail1 = expdata['00-10'].tail(1)
    tail2 = expdata['10-20'].tail(1)
    x_vals_1 = expdata['00-10']['x'] 
    x_vals_2 = expdata['10-20']['x']
    last_data_point_1 = [tail1['x'].iloc[0], tail1['y'].iloc[0]]
    last_data_point_2 = [tail2['x'].iloc[0], tail2['y'].iloc[0]]

    ## Read in the prompt photons
    prompt_loc = ROOT_DIRECTORY+"v2/jetscape_data/prompt_photons/AuAu_200/"
    prompt_0_10_tmp, kfactor_1  = pfuncs.get_prompt_photons(prompt_loc, multiplicity['00-10'], last_data_point_1, sigma)
    prompt_10_20_tmp, kfactor_2 = pfuncs.get_prompt_photons(prompt_loc, multiplicity['10-20'], last_data_point_2, sigma, kfactor_1)
    print(kfactor_1, kfactor_2)
    prompt_0_10  = pfuncs.harmonize_x_axis(prompt_0_10_tmp, x_vals_1, 'prompt')
    prompt_10_20 = pfuncs.harmonize_x_axis(prompt_10_20_tmp, x_vals_2, 'prompt')
    ## Thermal Photons:
    #THERMAL_ROOT_DIR = "/Users/rmyazdi/Documents/research/hydro_plots/VISHNU/HydroProfiles/"
    # thermal_0_10   = pfuncs.get_thermal_photons('AuAu200', '00-10', THERMAL_ROOT_DIR)
    # thermal_10_20  = pfuncs.get_thermal_photons('AuAu200', '10-20', THERMAL_ROOT_DIR)
    # thermal_0_10  = pfuncs.harmonize_x_axis(thermal_0_10, x_vals_1, channel='thermal')
    # thermal_10_20 = pfuncs.harmonize_x_axis(thermal_10_20, x_vals_2, channel='thermal')
    KOMPOST_ROOT_DIR = ROOT_DIRECTORY + 'v2/other_data/JF_MultiMessenger/'
    thermal_0_10 = pfuncs.get_JF_photons('AuAu200', '00-10', 'thermal', KOMPOST_ROOT_DIR)
    thermal_0_10 = pfuncs.harmonize_x_axis(thermal_0_10, x_vals_1, channel='thermal')
    thermal_10_20 = pfuncs.get_JF_photons('AuAu200', '10-20', 'thermal', KOMPOST_ROOT_DIR)
    thermal_10_20 = pfuncs.harmonize_x_axis(thermal_10_20, x_vals_2, channel='thermal')

    ## Pre-Equilibrium
    KOMPOST_ROOT_DIR = ROOT_DIRECTORY + 'v2/other_data/JF_MultiMessenger/'
    kompost_0_10 = pfuncs.get_JF_photons('AuAu200', '00-10', 'preEq', KOMPOST_ROOT_DIR)
    kompost_0_10 = pfuncs.harmonize_x_axis(kompost_0_10, x_vals_1, channel='kompost')
    kompost_10_20 = pfuncs.get_JF_photons('AuAu200', '10-20', 'preEq', KOMPOST_ROOT_DIR)
    kompost_10_20 = pfuncs.harmonize_x_axis(kompost_10_20, x_vals_2, channel='kompost')
    
    ## Jet Medium Photons:
    JMEDIUM_ROOT_DIR = ROOT_DIRECTORY + 'v2/jetscape_data/'
    martini_0_10   = pfuncs.get_jet_medium_photons('martini', '00-10', sigma, '200', 'AuAu200', multiplicity, JMEDIUM_ROOT_DIR)
    martini_10_20  = pfuncs.get_jet_medium_photons('martini', '10-20', sigma, '200', 'AuAu200', multiplicity, JMEDIUM_ROOT_DIR)
    martini_0_10  = pfuncs.harmonize_x_axis(martini_0_10, x_vals_1, channel='jet-medium')
    martini_10_20 = pfuncs.harmonize_x_axis(martini_10_20, x_vals_2, channel='jet-medium')

    cujet_0_10   = pfuncs.get_jet_medium_photons('cujet', '00-10', sigma, '200', 'AuAu200', multiplicity, JMEDIUM_ROOT_DIR)
    cujet_10_20  = pfuncs.get_jet_medium_photons('cujet', '10-20', sigma, '200', 'AuAu200', multiplicity, JMEDIUM_ROOT_DIR)
    cujet_0_10  = pfuncs.harmonize_x_axis(cujet_0_10, x_vals_1, channel='jet-medium')
    cujet_10_20 = pfuncs.harmonize_x_axis(cujet_10_20, x_vals_2, channel='jet-medium')

    ## Construct totals:
    totals = {}

    tot_no_jmed_y   = prompt_0_10['prompt'] + kompost_0_10['kompost'] + thermal_0_10['thermal'] 
    dtot_no_jmed_y  = prompt_0_10['dprompt']
    tot_martini_y   = tot_no_jmed_y + martini_0_10['total']
    dtot_martini_y  = np.sqrt(dtot_no_jmed_y** 2 + martini_0_10['dtotal']**2)
    tot_cujet_y     = tot_no_jmed_y + cujet_0_10['total']
    dtot_cujet_y    = np.sqrt(dtot_no_jmed_y** 2 + cujet_0_10['dtotal']**2)
    totals['00-10'] = pd.DataFrame({'pT':x_vals_1, 'No Jet-Med.':tot_no_jmed_y, 'err No Jet-Med.':dtot_no_jmed_y,
                                    'MATTER+MARTINI':tot_martini_y, 'err MATTER+MARTINI':dtot_martini_y,
                                    'MATTER+CUJET':tot_cujet_y, 'err MATTER+CUJET':dtot_cujet_y})

    tot_no_jmed_y   = prompt_10_20['prompt'] + kompost_10_20['kompost'] + thermal_10_20['thermal']
    dtot_no_jmed_y  = prompt_10_20['dprompt']
    tot_martini_y   = tot_no_jmed_y + martini_10_20['total']
    dtot_martini_y  = np.sqrt(dtot_no_jmed_y** 2 + martini_10_20['dtotal']**2)
    tot_cujet_y     = tot_no_jmed_y + cujet_10_20['total']
    dtot_cujet_y    = np.sqrt(dtot_no_jmed_y** 2 + cujet_10_20['dtotal']**2)
    totals['10-20'] = pd.DataFrame({'pT':x_vals_2, 'No Jet-Med.':tot_no_jmed_y, 'err No Jet-Med.':dtot_no_jmed_y,
                                    'MATTER+MARTINI':tot_martini_y, 'err MATTER+MARTINI':dtot_martini_y,
                                    'MATTER+CUJET':tot_cujet_y, 'err MATTER+CUJET':dtot_cujet_y})
    ## Construct ratios:
    ratios = {'00-10':{}, '10-20':{}}

    m_ratio_thermal_to_tot = thermal_0_10['thermal']/totals['00-10']['MATTER+MARTINI']
    m_ratio_kompost_to_tot = kompost_0_10['kompost']/totals['00-10']['MATTER+MARTINI']
    m_ratio_prompt_to_tot  = prompt_0_10['prompt']  /totals['00-10']['MATTER+MARTINI']
    m_ratio_jetmed_to_tot  = martini_0_10['total']  /totals['00-10']['MATTER+MARTINI']
    ratios['00-10']['martini'] = pd.DataFrame({'pT':x_vals_1,
                       'thermal':m_ratio_thermal_to_tot, 'jet-med':m_ratio_jetmed_to_tot,
                       'kompost':m_ratio_kompost_to_tot, 'prompt':m_ratio_prompt_to_tot})

    c_ratio_thermal_to_tot = thermal_0_10['thermal']/totals['00-10']['MATTER+CUJET']
    c_ratio_kompost_to_tot = kompost_0_10['kompost']/totals['00-10']['MATTER+CUJET']
    c_ratio_prompt_to_tot  = prompt_0_10['prompt']  /totals['00-10']['MATTER+CUJET']
    c_ratio_jetmed_to_tot  = cujet_0_10['total']    /totals['00-10']['MATTER+CUJET']
    ratios['00-10']['cujet'] = pd.DataFrame({'pT':x_vals_1,
                       'thermal':c_ratio_thermal_to_tot, 'jet-med':c_ratio_jetmed_to_tot,
                       'kompost':c_ratio_kompost_to_tot, 'prompt':c_ratio_prompt_to_tot})

    m_ratio_thermal_to_tot = thermal_10_20['thermal']/totals['10-20']['MATTER+MARTINI']
    m_ratio_kompost_to_tot = kompost_10_20['kompost']/totals['10-20']['MATTER+MARTINI']
    m_ratio_prompt_to_tot  = prompt_10_20['prompt']  /totals['10-20']['MATTER+MARTINI']
    m_ratio_jetmed_to_tot  = martini_10_20['total']  /totals['10-20']['MATTER+MARTINI']
    ratios['10-20']['martini'] = pd.DataFrame({'pT':x_vals_2,
                       'thermal':m_ratio_thermal_to_tot, 'jet-med':m_ratio_jetmed_to_tot,
                       'kompost':m_ratio_kompost_to_tot, 'prompt':m_ratio_prompt_to_tot})

    c_ratio_thermal_to_tot = thermal_10_20['thermal']/totals['10-20']['MATTER+CUJET']
    c_ratio_kompost_to_tot = kompost_10_20['kompost']/totals['10-20']['MATTER+CUJET']
    c_ratio_prompt_to_tot  = prompt_10_20['prompt']  /totals['10-20']['MATTER+CUJET']
    c_ratio_jetmed_to_tot  = cujet_10_20['total']    /totals['10-20']['MATTER+CUJET']
    ratios['10-20']['cujet'] = pd.DataFrame({'pT':x_vals_2,
                       'thermal':c_ratio_thermal_to_tot, 'jet-med':c_ratio_jetmed_to_tot,
                       'kompost':c_ratio_kompost_to_tot, 'prompt':c_ratio_prompt_to_tot})
    ## Plotting time
    labels =['No Jet-Med.','MATTER+MARTINI','MATTER+CUJET']
    alpha=0.3
    fig1, axes1 = plt.subplots(2, 2, height_ratios=(3,1), sharex=True, sharey='row')
    axes = axes1.flatten()

    for icent, centrality in enumerate(expdata):
        data = expdata[centrality]
        sims = totals[centrality]
        ax = axes[icent]
        ax.set_yscale('log')
        ax.text(0.05, 0.05, centrality + r'$\%$', transform=ax.transAxes)
        util.plot_expr_data_on_axis(ax, data, marker='s', color='black', face='black', s=20)
        for label in labels: 
            col = module_colors[label]
            y, dy = sims[label], sims[f'err {label}']
            ax.plot(sims['pT'], y, color=col, label=label)
            ax.fill_between(sims['pT'], y+dy, y-dy, color=col, alpha=alpha)
            deltax = 0.5*(data["xhigh"]-data["xlow"])
 
            ## form ratio to data and plot in the lower row
            
            ratio = y/data['y']
            err = dy**2 / y**2
            dratio_plus_stat = ratio * np.sqrt(err+ data['dratio_stat_plus'])
            dratio_minus_stat = ratio * np.sqrt(err + data['dratio_stat_minus'])
            dratio_plus_syst = ratio * data['dratio_syst_plus']
            dratio_minus_syst = ratio * data['dratio_syst_minus']
            errorboxes = [Rectangle((x-delx, y - yerrlow), width=2*delx, 
                                    height=abs(yerrlow)+abs(yerrhigh), zorder=5)
                            for x, delx, y, yerrlow, yerrhigh in
                            zip(data["x"], deltax, ratio, dratio_plus_syst, dratio_minus_syst)]
            pc = PatchCollection(errorboxes, facecolor=col, alpha=0.2)
            axes[icent+2].add_collection(pc)
            axes[icent+2].scatter(data["x"] , ratio, color=col, marker='p', s=10)
            axes[icent+2].errorbar(data["x"], ratio, xerr=deltax, yerr=[dratio_minus_stat, dratio_plus_stat], 
                                        fmt='none', 
                                        color=col)

    for ax in axes[2:]:
        ax.set_xlabel(r'$p_T$ (GeV)')
        ax.axhline(y=1., linestyle='dotted', color='black')

    handles = [Line2D([],[],color=c, label=l) for l, c in module_colors.items() if l != 'MATTER']
    handles.append(Line2D([],[],linewidth=None, color='black', marker='s', markersize=5, label=r'PHENIX (2012) $|\eta|<0.35$'))
    axes[1].legend(loc='upper right', handles=handles, fontsize=25)
    axes[2].set_ylabel('Theory over'+'\n data')
    axes[0].set_ylabel(r'$\frac{1}{2\pi\;p_T}\frac{dN}{dp_T\;d\eta}$ (GeV)${}^{-2}$')

    ## Now plot channel ratios to total:
    fig2, axes2 = plt.subplots(2, 2, sharex=True, sharey=True)
    for irow, centrality in enumerate(ratios):
        specs = ratios[centrality]
        for icol, eloss in enumerate(specs):
            spec = specs[eloss]
            ax = axes2[irow][icol]
            if irow == 0:
                ax.set_title(eloss.upper())
            for channel in ['thermal','prompt','kompost']:
                ax.plot(spec['pT'], spec[channel], color=channel_colors[channel])
            ax.plot(spec['pT'], spec['jet-med'], color=module_colors['MATTER+'+eloss.upper()])
            ax.text(0.05,0.8, centrality+r'$\%$', transform=ax.transAxes)
    for ax in axes2[1,:]:
        ax.set_xlabel(r'$p_T$ (GeV)')
    for i in [0,1]:
        axes2[i][0].set_ylabel(r'Channel to '+'\n'+ 'Total')
    handles = [Line2D([],[],color=channel_colors['prompt'], label='Prompt'),
               Line2D([],[],color=channel_colors['kompost'],label='Pre. Eq.'),
               Line2D([],[],color=channel_colors['thermal'],label='Thermal'),
               Line2D([],[],color=module_colors['MATTER+MARTINI'], label='MATTER+MARTINI'),
               Line2D([],[],color=module_colors['MATTER+CUJET'],label='MATTER+CUJET')]
    axes2[0][1].legend(loc='upper left', handles=handles, bbox_to_anchor=(1.05,0.9))
    plt.show()    
