from pandas import read_csv, DataFrame

def get_charged_spec(fname, pcut_low, pcut_high):
    tmp = read_csv(fname, comment='#')
    tmp['pT']  = 0.5 * (tmp['pTmin'] + tmp['pTmax'])
    tmp['dpT'] = tmp['pTmax'] - tmp['pTmin']
    spec_aft_cut = tmp[tmp['pT'].between(pcut_low,pcut_high)]
    return spec_aft_cut


def get_jet_spec(fname, pcut_low, pcut_high):
    tmp = read_csv(fname, comment='#')
    tmp['pT'] = 0.5 * (tmp['pTmin'] + tmp['pTmax'])
    tmp['dpT'] = tmp['pTmax'] - tmp['pTmin']
    result = tmp[tmp['pT'].between(pcut_low, pcut_high)]
    return result

def get_jet_shape(fname, rmax):
    calcs = read_csv(fname,comment='#')
    njet = 1
    with open(fname, 'r') as f:
        line = f.readline()
        line = line.split(' ')[-1]
        njet = float(line)

    dat = calcs[calcs['rmax'] < rmax]
    delta_r = dat['rmax'] - dat['rmin']
    r = 0.5*(dat['rmax'] + dat['rmin'])
    rho  = dat['N']  /(njet)
    drho = dat['dN'] /(njet)
    norm = sum(rho.to_list())#
    rho_normed = rho  /(delta_r * norm)
    drho_normed = drho/(delta_r * norm)
    r = r.to_list()
    dr = delta_r.to_list()
    rho_normed=rho_normed.to_list()
    drho_normed=drho_normed.to_list()
    del rho
    del drho
    del dat
    return DataFrame({'r':r,'dr':dr,'rho':rho_normed,'drho':drho_normed})
