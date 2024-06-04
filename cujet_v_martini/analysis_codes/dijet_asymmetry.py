import sys
import os
from Helpers import Histogram, get_xsection, Jet
import observables as obs

"""
    Do dijets only for the inclusive jets
"""
def collect_dijet_pTHatbin(fname, props, pt1vals):
    nevt = 0
    num_histograms = len(pt1vals)-1
    xmin = props['xj bins'][0]
    hists = [Histogram(props['xj bins']) for i in range(num_histograms)]
    ndijet = [0 for i in range(num_histograms)]
    jets_in_evt = []
    pt2min = props['xj bins'][0]*pt1vals[0]
    for line in open(fname, 'r'):
        line = line.split(',')
        len_line = len(line)
        if len_line == 1:
            nevt += 1 
        elif len_line == 5:
            currjet = Jet(line)
            if abs(currjet.eta) < props['eta cut']:
                jets_in_evt.append(currjet)
            if len(jets_in_evt) >= 2:
                jet_1 = jets_in_evt[0]
                jet_2 = jets_in_evt[1]
                if jet_1.calc_del_phi(jet_2) > props['phi min']:
                    for i in range(num_histograms):
                        if pt1vals[i] <= jet_1.pT < pt1vals[i+1] and pt1vals[i]*xmin <= jet_2.pT:
                            hists[i].fill(jet_2.pT/jet_1.pT)
                            ndijet[i] += 1
                jets_in_evt.clear()
                    
        else:
            continue
    return nevt,ndijet,hists

def average_pTHat(dirname, nruns, jet_radius, props, pt1vals):
    num_histograms = len(pt1vals)-1
    xmin           = props['xj bins'][0]
    dijets         = [Histogram(props['xj bins']) for i in range(num_histograms)]
    num_dijets     = [0 for i in range(num_histograms)]
    num_evts = 0
    for irun in range(nruns):
        curr_dir   = dirname.format(irun)
        temp_xsec  = get_xsection(curr_dir)
        fname_ncut = curr_dir + f"FinalStateJets_{jet_radius:0.6f}_ncut.dat" 
        tmp_nevt, tmp_num_dijet, tmp_dijets = collect_dijet_pTHatbin(fname_ncut, props, pt1vals)
        
        num_evts += tmp_nevt
        for i in range(num_histograms):
            tmp_dijets[i].multiply_by(temp_xsec)
            dijets[i] += tmp_dijets[i]
            num_dijets[i] += tmp_num_dijet[i]
    ## TODO correct nevt for dijet 
    #scale = 1./nruns
    #num_dijets = [e*scale for e in num_dijets]
    #for hist in dijets:
    #    hist.multiply_by(scale)
    return num_evts, num_dijets, dijets

def build_dijet_spec(dirname, pthatbins, nruns, jet_radius, props, pt1vals):
    num_histograms = len(pt1vals)-1
    num_pThat_bins = len(pthatbins)-1
    xmin           = props['xj bins'][0]
    dijets         = [Histogram(props['xj bins']) for i in range(num_histograms)]
    num_dijets     = [0 for i in range(num_histograms)]
    nevts = 0
    for ipT in range(num_pThat_bins):
        ptHatmin, ptHatmax = pthatbins[ipT], pthatbins[ipT+1]
        pthatdir = dirname + f"{ptHatmin}_{ptHatmax}_" + "{0:d}/"
        tmp_nevt, tmp_ndijets, tmp_dijets = average_pTHat(pthatdir, nruns, jet_radius, props, pt1vals)
        nevts += tmp_nevt
        for i in range(num_histograms):
            dijets[i] += tmp_dijets[i]
            num_dijets[i] += tmp_ndijets[i]

    return nevts, num_dijets, dijets

def save_to_file(fname, props, data, nevts):
    """
        pt1vals, spectrum, ndijet
        data = {(ptmin, ptmax): (num dijets, dijet spec)}
    """
    xbins = props['xj bins']
    eta_cut = props['eta cut']
    phi_cut = props['phi min']
    with open(fname, 'w') as f:
        f.write(f"#eta_cut {eta_cut:0.2f} del_phi {phi_cut:0.5f} nevts {nevts}\n")
        header = 'xmin,xmax'
        for key in data.keys():
            header += f',{key[0]:d}-{key[1]:d},err {key[0]:d}-{key[1]:d}'
        header += '\n'
        f.write(header)
        ## write the number of dijets as a the first row
        f.write("-1,-1")
        for key in data.keys():
            f.write(f",{data[key][0]:0.5e},0")
        f.write("\n")
        for ix in range(len(xbins)-1):
            xmin, xmax = xbins[ix], xbins[ix+1]
            f.write(f"{xmin:0.2f},{xmax:0.2f}")
            for key in data.keys():
                val, err = data[key][1][ix]
                f.write(f",{val:0.5e}-{err:0.5e}")
            f.write("\n")
    print(f"Done writing {fname}")

if __name__=='__main__': 
    ## get the system
    system = sys.argv[1]
    eloss = ''
    centrality = ''
    maxT = ''

    location = '/home/rmyazdi/scratch/JETSCAPE_local/jetscape_exec/gamma_proj/'
    if 'pp' not in system:
        eloss      = sys.argv[2]
        centrality = sys.argv[3]
        location += f'HIC/{eloss}/{system}_{centrality}/'
    else:
        maxT = sys.argv[2]
        location += f'pp/results/maxT_{maxT}/{system}/'

    index_of_system = -1
    if '200' in system:
        index_of_system = 0
    elif '2760' in system:
        index_of_system = 1
    else:
        index_of_system = 2

    tag = 'new'
    if system == 'PbPb2760' and centrality in ['00-05','20-30','30-40']:
        if eloss != 'cujet':
            tag = 'old' 
    
    hard_prop = obs.dijet_hard_prop
    hard_pT   = obs.dijet_hard_pT

    miniJet_prop = obs.dijet_minijet_prop
    miniJet_pT   = obs.dijet_minijet_pT

    pThats = obs.pT_hat_bins[index_of_system]
    nruns = obs.nruns[tag]

    save_name = '/home/rmyazdi/scratch/JETSCAPE_local/jetscape_exec/gamma_proj/results/'
    if maxT != '':
        save_name+= f'maxT_{maxT}/'
    if eloss != '':
        save_name += f'{eloss}/'
    save_name += f'{system}_'
    if centrality != '':
        save_name += f'{centrality}_'

    jet_radii = [0.4]#[0.2, 0.3, 0.4]#, 0.6, 0.8]
    print("PROCESS DIJETS")
    for R in jet_radii:
        print(f"Jet Radius: {R}")
        results = {}
        ## for the two LHC systems:
        if index_of_system != 0:
            (nevts, ndijet, hist) = build_dijet_spec(location, pThats, nruns, R, hard_prop, hard_pT)
            for ipt in range(len(hard_pT)-1):
                ptmin, ptmax = hard_pT[ipt], hard_pT[ipt+1]
                results[(ptmin, ptmax)] = (ndijet[ipt], hist[ipt])
        
            dijet_spec_fname = save_name + f'dijet_spec_hard_jets_radius_{R}.csv'
            save_to_file(dijet_spec_fname, hard_prop, results, nevts)
            
        del results
        results = {}
        (nevts, ndijet, hist) = build_dijet_spec(location, pThats, nruns, R, miniJet_prop, miniJet_pT)
        for ipt in range(len(miniJet_pT)-1):
            ptmin, ptmax = miniJet_pT[ipt], miniJet_pT[ipt+1]
            results[(ptmin, ptmax)] = (ndijet[ipt], hist[ipt])
         
        dijet_spec_fname = save_name + f'dijet_spec_mini_jets_radius_{R}.csv'
        save_to_file(dijet_spec_fname, miniJet_prop, results, nevts)
