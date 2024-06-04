import numpy as np
import pandas as pd

"""
    One single plot for jet shape ratio and 
    jet fragmentation function as a function of 
    z \equiv pt(jet).pt(had)/pt(jet)**2
    No need for more detail. 
"""
zmin, zmax = 0.005, 2
## JET SHAPE STUFF
## My new MARTINI runs with and without jet formation time
def read_martini_jet_shape(fname):
    tmp = pd.read_csv(fname,comment='#')
    njet = 1
    with open(fname, 'r') as f:
        line = f.readline()
        line = line.split(' ')[-1]
        njet = float(line)
    dat     = tmp[tmp['rmax'] < 0.31]
    delta_r = tmp['rmax'] - tmp['rmin']
    r       = 0.5*(tmp['rmax'] + tmp['rmin'])
    rho  = tmp['N']  /(njet)
    drho = tmp['dN'] /(njet)
    norm = sum(rho.to_list())#
    rho_normed  = rho /(delta_r * norm)
    drho_normed = drho/(delta_r * norm)
    return rho_normed, drho_normed, r, delta_r

def get_shape_with_formation_time(path_pp, path_AA):
    ## Get the Proton-Proton Baseline for the case with formation time:
    path = path_pp + "jet_shape.csv" # ../../calcs/pp_2p76/
    calcs = pd.read_csv(path,comment='#')
    njet = 1
    with open(path, 'r') as f:
        line = f.readline()
        line = line.split(' ')[-1]
        njet = float(line)

    dat = calcs[calcs['rmax'] < 0.31]
    delta_r = dat['rmax'] - dat['rmin']
    r = 0.5*(dat['rmax'] + dat['rmin'])
    rho  = dat['N']  /(njet)
    drho = dat['dN'] /(njet)
    norm = sum(rho.to_list())#
    rho_normed = rho  /(delta_r * norm)
    drho_normed = drho/(delta_r * norm)

    pp_data = pd.DataFrame({'r':r, 'dr':delta_r, 
                            'rho_normed':rho_normed, 
                            'drho_normed':drho_normed})

    ## Compute the jet shape ratio for each rate set and the two 
    ## centralities of interest:
    fname = path_AA + "jet_shape.csv" # ../../calcs/final_runs/rset_1/2p76/{c}
    rho_normed, drho_normed, r, delta_r = read_martini_jet_shape(fname.format(c='0_5'))
    scale = sum(rho_normed)
    rho_normed = rho_normed/(delta_r*scale)
    drho_normed = drho_normed/(delta_r*scale)
    ratio = rho_normed/pp_data['rho_normed']
    dr = ratio*np.sqrt((drho_normed/rho_normed)**2 + (pp_data['drho_normed']/pp_data['rho_normed'])**2)
    return pd.DataFrame({'r':r, 'delta_r':delta_r, 'ratio':ratio, 'dratio':dr})

def get_shape_without_formation_time(path_pp, path_AA):
    ## Get the Proton-Proton Baseline for the case with formation time:
    path = path_pp + "jet_shape.csv" 
    calcs = pd.read_csv(path,comment='#')
    njet = 1
    with open(path, 'r') as f:
        line = f.readline()
        line = line.split(' ')[-1]
        njet = float(line)

    dat = calcs[calcs['rmax'] < 0.31]
    delta_r = dat['rmax'] - dat['rmin']
    r = 0.5*(dat['rmax'] + dat['rmin'])
    rho  = dat['N']  /(njet)
    drho = dat['dN'] /(njet)
    norm = sum(rho.to_list())#
    rho_normed = rho  /(delta_r * norm)
    drho_normed = drho/(delta_r * norm)

    pp_data = pd.DataFrame({'r':r, 'dr':delta_r, 
                            'rho_normed':rho_normed, 
                            'drho_normed':drho_normed})

    ## Compute the jet shape ratio for each rate set and the two 
    ## centralities of interest:
    # fname = "../../../v2/production_v2/martini_results/final_PbPb_2p76"
    # fname += "/rset_1/cent_{c}/
    fname = path_AA+"jet_shape.csv"
    rho_normed, drho_normed, r, delta_r = read_martini_jet_shape(fname.format(c='0_5'))
    scale = sum(rho_normed)
    rho_normed = rho_normed/(delta_r*scale)
    drho_normed = drho_normed/(delta_r*scale)
    ratio = rho_normed/pp_data['rho_normed']
    dr = ratio*np.sqrt((drho_normed/rho_normed)**2 + (pp_data['drho_normed']/pp_data['rho_normed'])**2)
    return pd.DataFrame({'r':r, 'delta_r':delta_r, 'ratio':ratio, 'dratio':dr})

## JETSCAPE RUNS FOR JET SHAPE:
def jetscape_pp_jet_shape(fname):
    ## Read in the PP jet shape and massage it:
    #path_pp = '../../jetscape_data/max_time/maxT_200_highstat/pp_2760_jet_shape_LHC_cuts_jet_rad_0.3_0.30_2.00.csv'
    tmp = pd.read_csv(fname,comment='#')
    nevt_wcut, njet_wcut = 1, 1
    with open(fname, 'r') as f:
        line = f.readline()
        line = line.split()
        (nevt_wcut, njet_wcut) = float(line[-3]), float(line[-1])
    tmp = tmp[tmp['rmax'] < 0.31]
    delta_r = tmp['rmax'] - tmp['rmin']
    r = 0.5*(tmp['rmax'] + tmp['rmin'])
    rho  = tmp['wcut']  /(nevt_wcut * njet_wcut )
    drho = tmp['dwcut'] /(nevt_wcut * njet_wcut )
    norm = sum(rho.to_list())#
    rho_normed = rho  /(delta_r * norm)
    drho_normed = drho/(delta_r * norm)
    return pd.DataFrame({'r':r, 'dr':delta_r, 'rho':rho_normed, 'drho':drho_normed})

def jetscape_aa_jet_shape(fname):

    tmp = pd.read_csv(fname, comment='#')
    nevt_wcut, njet_wcut = 1, 1
    with open(fname, 'r') as f:
        line = f.readline()
        line = line.split()
        (nevt_wcut, njet_wcut) = float(line[-3]), float(line[-1])
    tmp = tmp[tmp['rmax'] < 0.31]
    delta_r = tmp['rmax'] - tmp['rmin']
    r = 0.5*(tmp['rmax'] + tmp['rmin'])
    rho  = tmp['wcut']  /(nevt_wcut * njet_wcut )
    drho = tmp['dwcut'] /(nevt_wcut * njet_wcut )
    norm = sum(rho.to_list())
    rho_normed = rho  /(delta_r * norm)
    drho_normed = drho/(delta_r * norm)
    return pd.DataFrame({'r':r, 'dr':delta_r, 'rho':rho_normed, 'drho':drho_normed})

def process_jetscape_jet_shape_results(pp_fname, aa_fname):
    pp = jetscape_pp_jet_shape(pp_fname)
    aa = jetscape_aa_jet_shape(aa_fname)

    ratio = aa['rho']/pp['rho']
    dratio = ratio * np.sqrt( (aa['drho']/aa['rho'])**2 + (pp['drho']/pp['rho'])**2 )
    return pd.DataFrame({'r':pp['r'], 'dr':pp['dr'], 'ratio':ratio, 'dratio':dratio})


"""
    JET FRAGMENTATION FUNCTION RATIOS
"""
def get_FFz_martini_results_wformtime(fname):
    #fname_z  = f'../../calcs/final_runs/rset_1/2p76/0_5/jet_FF_z.csv'
    tmp_DZ   = pd.read_csv(fname,comment='#')
    tmp_DZ   = tmp_DZ[tmp_DZ['zmax'] < zmax]
    tmp_DZ   = tmp_DZ[tmp_DZ['zmin'] > zmin]
    njets_pp, nevts_pp = 0 , 0
    with open(fname,'r') as f:
        line = f.readline()
        line = line.split(' ')
        njets_pp = float(line[1])
        # nevts_pp = float(line[-1])
    scale = njets_pp
    dz = tmp_DZ['zmax'] - tmp_DZ['zmin']
    z = 0.5*(tmp_DZ['zmax'] + tmp_DZ['zmin'])
    tmpyz  = tmp_DZ['N'] /(dz*scale)
    dtmpyz = tmp_DZ['dN']/(dz*scale)
    result_z = pd.DataFrame({'z':z,
                           'dz':dz,
                           'FF' :tmpyz,
                           'dFF':dtmpyz})
    return result_z

def get_FFz_martini_results_nformtime(fname):
    #fname_z  = f'../../../v2/production_v2/martini_results/final_PbPb_2p76/rset_{rateset}/cent_{cent}/jet_FF_z.csv'

    tmp_DZ   = pd.read_csv(fname, comment='#')
    tmp_DZ   = tmp_DZ[tmp_DZ['zmax'] < zmax]
    tmp_DZ   = tmp_DZ[tmp_DZ['zmin'] > zmin]
    njets_pp, nevts_pp = 0 , 0
    with open(fname, 'r') as f:
        line = f.readline()
        line = line.split(' ')
        njets_pp = float(line[1])
    
    scale = njets_pp
    dz = tmp_DZ['zmax'] - tmp_DZ['zmin']
    z = 0.5*(tmp_DZ['zmax'] + tmp_DZ['zmin'])
    tmpyz  = tmp_DZ['N'] /(dz*scale)
    dtmpyz = tmp_DZ['dN']/(dz*scale)
    result_z = pd.DataFrame({'z':z,
                           'dz':dz,
                           'FF' :tmpyz,
                           'dFF':dtmpyz})
    return result_z

def get_FFz_pp_wformtime(fname):
    ## Read in the pp data:
    pp_DZ      = pd.read_csv(fname,comment='#')

    pp_DZ = pp_DZ[pp_DZ['zmax'] < zmax]
    pp_DZ = pp_DZ[pp_DZ['zmin'] > zmin]
    njets_pp, nevts_pp = 0 , 0
    with open(fname,'r') as f:
        line = f.readline()
        line = line.split(' ')
        njets_pp = float(line[1])
        # nevts_pp = float(line[-1])

    scale = njets_pp
    dz = pp_DZ['zmax'] - pp_DZ['zmin']
    z = 0.5*(pp_DZ['zmax'] + pp_DZ['zmin'])
    y  = pp_DZ['N'] /(dz*scale)
    dy = pp_DZ['dN']/(dz*scale)
    pp_dz = pd.DataFrame({'z':z,'dz':dz,'FF':y,'dFF':dy})
    return pp_dz

def get_FFz_pp_nformtime(fname):
    ## Read in the pp data:
    pp_DZ      = pd.read_csv(fname,comment='#')

    pp_DZ = pp_DZ[pp_DZ['zmax'] < zmax]
    pp_DZ = pp_DZ[pp_DZ['zmin'] > zmin]
    njets_pp, nevts_pp = 0 , 0
    with open(fname, 'r') as f:
        line = f.readline()
        line = line.split(' ')
        njets_pp = float(line[1])
        nevts_pp = float(line[-1])

    scale = njets_pp
    dz = pp_DZ['zmax'] - pp_DZ['zmin']
    z = 0.5*(pp_DZ['zmax'] + pp_DZ['zmin'])
    y  = pp_DZ['N'] /(dz*scale)
    dy = pp_DZ['dN']/(dz*scale)
    pp_dz = pd.DataFrame({'z':z,'dz':dz,'FF':y,'dFF':dy})
    return pp_dz

def get_FFz_ratio_martini_with_form_ratio(pp_fname, aa_fname):
    AA = get_FFz_martini_results_wformtime(aa_fname)
    pp = get_FFz_pp_wformtime(pp_fname)

    z, dz = AA['z'], AA['dz']
    r = AA['FF']/pp['FF']
    dr = r * np.sqrt((AA['dFF']/AA['FF'])**2 + (pp['dFF']/pp['FF'])**2)
    return pd.DataFrame({'z':z, 'dz':dz, 'ratio':r, 'dratio':dr})

def get_FFz_ratio_martini_without_form_ratio(pp_fname, aa_fname):
    AA = get_FFz_martini_results_nformtime(aa_fname)
    pp = get_FFz_pp_nformtime(pp_fname)

    z, dz = AA['z'], AA['dz']
    r = AA['FF']/pp['FF']
    dr = r * np.sqrt((AA['dFF']/AA['FF'])**2 + (pp['dFF']/pp['FF'])**2)
    return pd.DataFrame({'z':z, 'dz':dz, 'ratio':r, 'dratio':dr})
    


## JETSCAPE RUNS FOR JET FRAGMENTATION FUNCTION
def get_FFz_pp_jetscape(fname):
    data = pd.read_csv(fname, comment='#')
    data = data[data['zmax'] < zmax]
    data = data[data['zmin'] > zmin]
    pp_jetspec = pd.read_csv(fname.replace('FF_z','jet_spec_FF') ,comment='#')
    njets_pp =  sum(pp_jetspec['Ncut'].tolist())

    z = 0.5 * (data['zmax'] + data['zmin'])
    dz = data['zmax'] - data['zmin']
    y = data['Ncut']/njets_pp
    dy = data['dNcut']/njets_pp

    return pd.DataFrame({'z':z, 'dz':dz, 'FF':y, 'dFF':dy})

def get_FFz_AA_jetscape(fname):
    data = pd.read_csv(fname, comment='#')
    print(data.columns)
    data = data[data['zmax'] < zmax]
    data = data[data['zmin'] > zmin]
    pp_jetspec = pd.read_csv(fname.replace('FF_z','jet_spec_FF') ,comment='#')
    njets_pp =  sum(pp_jetspec['Ncut'].tolist())

    z = 0.5 * (data['zmax'] + data['zmin'])
    dz = data['zmax'] - data['zmin']
    y = data['Ncut']/njets_pp
    dy = data['dNcut']/njets_pp

    return pd.DataFrame({'z':z, 'dz':dz, 'FF':y, 'dFF':dy})

def process_jet_FF_ratio_jetscape(fname_pp, fname_AA):
    pp = get_FFz_pp_jetscape(fname_pp)
    AA = get_FFz_AA_jetscape(fname_AA)

    z, dz = AA['z'], AA['dz']
    r = AA['FF']/pp['FF']
    dr = r * np.sqrt((AA['dFF']/AA['FF'])**2 + (pp['dFF']/pp['FF'])**2)
    return pd.DataFrame({'z':z, 'dz':dz, 'ratio':r, 'dratio':dr})