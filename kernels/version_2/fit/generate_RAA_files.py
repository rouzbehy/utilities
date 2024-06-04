import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.pyplot import cm
from matplotlib.lines import Line2D
import pandas as pd
import numpy as np
from matplotlib.patches import Rectangle
from matplotlib.collections import PatchCollection
from matplotlib import ticker as mticker
import dictionaries as ddicts
import jetDicts
from scipy.interpolate import interp1d
import info
import util

plt.rcParams.update(info.my_rcParams)


"""
This module just generates the RAA files
"""
kappas = info.kappas
all_kappas = [val for (k,val) in kappas.items()]
kappas_r, kappas_e = [e[0] for e in all_kappas],[e[1] for e in all_kappas]
kappa_sets = np.arange(0, 50)
pt_chgd_low = util.pt_chgd_low
pt_chgd_high = util.pt_chgd_high
pt_jet_low = util.pt_jet_low
pt_jet_high = util.pt_jet_high
jet_cone_radii = util.jet_cone_radii

if __name__=='__main__':
    ## Charged Hadron RAA data
    ALICE_HADS_00_05 = pd.read_csv("../../../exp_data/charged/ALICE/Table_16.csv", comment='#').rename(columns=ddicts.colnames_raa_alice)
    ALICE_HADS_00_05 = ALICE_HADS_00_05[ALICE_HADS_00_05['x'].between(pt_chgd_low,pt_chgd_high)]
    CMS_HADS_00_05   = pd.read_csv("../../../exp_data/charged/CMS/Table5_1.csv"  , comment='#').rename(columns=ddicts.colnames_raa_cms)
    CMS_HADS_00_05 = CMS_HADS_00_05[CMS_HADS_00_05['x'].between(pt_chgd_low,pt_chgd_high)]
    ATLAS_HADS_00_05 = pd.read_csv("../../../exp_data/charged/ATLAS/Table33.csv" , comment='#').rename(columns=ddicts.colnames_raa_atlas)
    ATLAS_HADS_00_05 = ATLAS_HADS_00_05[ATLAS_HADS_00_05['x'].between(pt_chgd_low,pt_chgd_high)]
    ## Jet RAA data
    CMS_JETS_00_05 = {f'0p{i}': pd.read_csv(f"../../../exp_data/jets/CMS/Table_JETRAA_0p{i}.csv",comment='#').rename(columns=jetDicts.colnames_CMS_RAA) for i in [2,3,4]}
    ## Jet Shape Ratio
    fname = '../../../exp_data/jet_shape/CMS_PbPb_2760_Jet_ShapeRatio_00-10_pT=100-INF_R=0p3_eta=0p3-2p0.dat'
    tmp                = np.loadtxt(fname,comments='#',unpack=True, delimiter='\t')
    CMS_JET_SHAPE      = pd.DataFrame({'x':tmp[0],'y':tmp[1],'dx':tmp[2],'dy':tmp[3],'xlow':tmp[0]-tmp[2]/2,'xhigh':tmp[0]+tmp[2]/2})

    ## Proton-Proton Baselines
    pp_chgd = util.get_charged_hadron_spec(None,None,'pp',pt_chgd_low, pt_chgd_high)
    pp_jets = util.get_jet_spec(None, None, 'pp', pt_jet_low, pt_jet_high)
    pp_jshp = util.get_jet_shape(None, None, 'pp')

    ## Pb-Pb calculations
    raa_calc_chgd = {}
    raa_calc_jet = {}
    raa_calc_shape = {}
    for rset in util.rset_shorts:
        raa_calc_chgd[rset] = {}
        raa_calc_jet[rset] = {}
        raa_calc_shape[rset] = {}
        for kset in kappa_sets:
            ## Charged hadrons:
            tmp = util.get_charged_hadron_spec(rset, kset, 'pbpb', pt_chgd_low, pt_chgd_high)
            tmp_raa = tmp['N']/pp_chgd['N']
            tmp_draa = tmp_raa*np.sqrt((pp_chgd['dN']/pp_chgd['N'])**2 + (tmp['dN']/tmp['N'])**2)
            raa_calc_chgd[rset][kset] = [pp_chgd['pT'], tmp_raa, tmp_draa]
            ## Jet Shape:
            tmp = util.get_jet_shape(rset, kset, 'pbpb')
            tmp_raa = tmp['rho_normed']/pp_jshp['rho_normed']
            tmp_draa = tmp_raa*np.sqrt((pp_jshp['drho_normed']/pp_jshp['rho_normed'])**2 + (tmp['drho_normed']/tmp['rho_normed'])**2)
            raa_calc_shape[rset][kset]  = [pp_jshp['r'], tmp_raa, tmp_draa]

            ## Jets (gotta loop over the cone radii)

            tmp = util.get_jet_spec(rset, kset, 'pbpb', pt_jet_low, pt_jet_high)
            raa_calc_jet[rset][kset] = {}
            for r in jet_cone_radii:
                raa  = tmp[f'N{r}']/pp_jets[f'N{r}']
                a, da = f'N{r}', f'dN{r}' 
                draa = raa*np.sqrt((tmp[da]/tmp[a])**2 + (pp_jets[da]/pp_jets[a])**2)
                raa_calc_jet[rset][kset][r] = (tmp['pT'],raa, draa)

    prefix='/Users/rmyazdi/Documents/research/KERNELS_NLO_NP_PART2/v3/calcs/RAA_for_fit/'

    jet_tmpl = '{:0.3f},{:0.6e},{:0.6e},{:0.6e},{:0.6e},{:0.6e},{:0.6e}'
    for rset in raa_calc_chgd:
        for kset in raa_calc_chgd[rset]:
            ## Start with Charged RAA
            with open(prefix+f'chgd_raa_rset_{rset}_kset_{kset}.csv', 'w') as f:
                f.write('pT,raa,draa\n')
                spec = raa_calc_chgd[rset][kset]
                for item in zip(*spec):
                    f.write(f"{item[0]:0.4f},{item[1]:0.6e},{item[2]:0.6e}\n")
            ## Jet Shape Ratio
            with open(prefix+f'jet_shape_ratio_rset_{rset}_kset_{kset}.csv', 'w') as f:
                f.write('r,Rrho,dRrho\n')
                spec = raa_calc_shape[rset][kset]
                for item in zip(*spec):
                    f.write(f"{item[0]:0.4f},{item[1]:0.6e},{item[2]:0.6e}\n")
            ## Jet RAA 
            with open(prefix+f"jet_raa_rset_{rset}_kset_{kset}.csv", 'w') as f:
                f.write("pT,r0p2,dr0p2,r0p3,dr0p3,r0p4,dr0p4\n")
                spec = raa_calc_jet[rset][kset]
                x,y1,dy1 = spec['0p2']
                _,y2,dy2 = spec['0p3']
                _,y3,dy3 = spec['0p4']
                for v in zip(x,y1,dy1,y2,dy2,y3,dy3):
                    #line = ','.join([str(item) for item in v])
                    line = jet_tmpl.format(*v)
                    f.write(line+"\n")
