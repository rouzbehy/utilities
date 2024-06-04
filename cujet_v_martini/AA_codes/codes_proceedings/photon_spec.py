#!/usr/bin/env python3
import numpy as np
import matplotlib as mpl
import matplotlib.ticker as plticker
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.collections import PatchCollection
from matplotlib.patches import Rectangle
from scipy.interpolate import interp1d
import util
import helper
import jetDicts as ddicts
from COLORS import module_colors, colors
my_rcParams={
    "text.usetex": True,
    "font.family": "Georgia",
    "font.size": 30,
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
    "axes.spines.top"  : False,
    "legend.frameon"   : False}
plt.rcParams.update(my_rcParams)


inel_Xsec = 62.03948 
elosses = {'cujet':'MATTER+CUJET', 'martini':'MATTER+MARTINI', 'none':'No Jet-Medium'}
label_fucks = {'cujet':'Spec. 3', 'martini':'Spec. 2', 'none':'Spec. 1'}
nice_channels = {'prompt':'Prompt','thermal':'Thermal','preEq':'Pre-Equilibrium'}
zorders = { 'cujet':0.5, 'none':0., 'martini':1.0 }

if __name__ == '__main__':
    expdata, pT_lower_lim, pT_upper_lim = helper.get_experimental_data()
    ratios, totals = helper.get_jetscape_data(expdata, pT_lower_lim, pT_upper_lim, xsec=inel_Xsec)

    fig, axes = plt.subplots(3,2,figsize=(16,9),sharey='row', sharex=True, height_ratios=(2,1,1))

    for icent, cent in enumerate(expdata.keys()):
        spec_ax = axes[0][icent]
        ratio_to_data_ax = axes[1][icent]
        ratio_to_total_ax = axes[2][icent]
        data = expdata[cent]
        util.plot_expr_data_on_axis(spec_ax, data, marker='o')
        spec_ax.text(0.1,0.1,f'{cent}\%', transform=spec_ax.transAxes)
        for eloss in elosses:
            pT, spec, dspec = totals[eloss][cent]
            col  = module_colors[elosses[eloss]]
            mark = ddicts.eloss_marker[eloss]
            zorder = zorders[eloss]
            spec_ax.plot(pT, spec, color=col, zorder=zorder)
            spec_ax.fill_between(pT, spec-dspec, spec+dspec, color=col, alpha=0.2, zorder=zorder)

            helper.take_ratio_to_data_and_plot(ratio_to_data_ax, data, [spec, dspec], col, mark, zorder)

            ## now plot the ratio of jetmed to total and othr to total
            x, ratio_jet_to_total, ratio_other_to_total = ratios[eloss][cent]
            ratio_to_total_ax.plot(x, ratio_jet_to_total, color=col, linestyle='solid')
            #ratio_to_total_ax.plot(x, ratio_other_to_total, color=col, linestyle='dotted')




    handles = [Line2D([],[],label=r'ALICE', marker='o',color='black',linestyle='none',markersize=10)]
    for eloss in elosses:
         label = elosses[eloss]
         handles.append(Line2D([],[],color=module_colors[label],label=label,marker=ddicts.eloss_marker[eloss],linestyle='none',markersize=10))
    axes[0][0].set_yscale('log')
    axes[0][0].set_ylabel(r'$\frac{1}{N_{\mathrm{evt}}\,2\pi\,p_T}\frac{\mathrm{d}N^{\gamma}}{\mathrm{d}p_T\mathrm{d}\eta}$, (GeV$^{-2}$)', fontsize=25)
    axes[0][0].set_ylim(top=1e0, bottom=1e-7)

    axes[1][0].set_ylabel(r'Total over'+'\n Data', fontsize=20 )
    axes[2][0].set_ylabel(r'Jet. Med. over'+'\n Total', fontsize=20)
    for ax in axes[1]:
        ax.axhline(1,color='black', linestyle='dashed')
    for ax in axes[2]:
        ax.set_xlabel(r'$p_T$ (GeV)', fontsize=30)
    other_handles = [Line2D([],[],color='black',linestyle='solid' ,label=r'$r_{\mathrm{j.med.}}$')]
                     #Line2D([],[],color='black',linestyle='dotted',label=r'$r_{\mathrm{other}}$')]
    #axes[2][1].legend(loc='upper left', handles=other_handles, fontsize=30, ncol=1)
    axes[0][1].legend(loc='upper right', handles=handles, fontsize=30)
    axes[0][0].text(0.25, 0.80, r'Pb-Pb @ $\sqrt{s}=2.76$ ATeV', transform=axes[0][0].transAxes)
    axes[0][0].text(0.40, 0.70, r'$|\eta|<0.8$',transform=axes[0][0].transAxes)
    # #fig.savefig("../../Plots/Photon_Plots/fig1_photon_spec_JF_prompts.pdf", dpi=200)
    plt.show()
