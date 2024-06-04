from matplotlib import colormaps as cmap

"""
    Define all the colors that I will use:
"""
tab_colors = cmap.get_cmap('tab10')
colors = tab_colors.colors

## LO/ NLO / NP Study:
#rate_set_colors = {'LO':colors[0], 'NLO':colors[1], 'NP':colors[2]}
## for the paper on LO, NLO, NP kernels,
rate_set_colors = {'LO': '#0077BB', 'NLO':'#EE7733', 'NP':'#009988'}
channel_colors = {'Thermal': colors[3], 'Prompt': colors[4], 'Pre-Eq.':colors[5]}



## 
module_colors = {'MATTER+MARTINI':colors[0], 
                 'CUJET+MATTER':colors[1], 
                 'MATTER':colors[2]}

maxt_colors = {'MATTER($\tau_{\mathrm{max}}=20$ fm/$c$)': colors[6],
               'MATTER($\tau_{\mathrm{max}}=200$ fm/$c$)': colors[7],
               'MATTER($\tau_{\mathrm{max}}=20$ fm/$c$)+MARTINI': colors[0],
               'MATTER($\tau_{\mathrm{max}}=20$ fm/$c$)+CUJET':colors[1],
               'MATTER($\tau_{\mathrm{max}}=200$ fm/$c$)+MARTINI':colors[8]}