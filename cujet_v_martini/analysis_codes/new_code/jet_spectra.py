import sys
import os
from Helpers import Histogram, get_xsection, Jet
import observables as obs

def collect_jets_pTHatbin(fname, jet_eta_cut, ptbins):
    """
        def collect_jets_pTHatbin(fname, jet_eta_cut, ptbins):
            read the file with _fname_ and populate 
            the `jet_hist` object with jets within
            the given `eta_cut` window
    """
    nevt = 0
    jet_hist = Histogram(ptbins)
    for line in open(fname,'r'):
        line = line.split(',')
        len_line = len(line)
        if len_line == 1:
            nevt += 1
            continue
        elif len_line == 5:## it's a jet
            currjet = Jet(line)
            if abs(currjet.eta) < jet_eta_cut:
                jet_hist.fill(currjet.pT)
        else:
            continue
    return nevt, jet_hist

def average_pTHat(dirname, nruns, eta_cut, ptbins, jet_radius):
    jet_hist_wcut = Histogram(ptbins)
    jet_hist_ncut = Histogram(ptbins)
    #jet_hist_chgd = Histogram(ptbins)

    nevt_ncut, nevt_wcut= 0,0
    for irun in range(nruns):
        curr_dir = dirname.format(irun)
        temp_xsec = get_xsection(curr_dir)

        fname_wcut = curr_dir + f"FinalStateJets_{jet_radius:0.6f}_wcut.dat" 
        fname_ncut = curr_dir + f"FinalStateJets_{jet_radius:0.6f}_ncut.dat" 
        #fname_chgd = curr_dir + f"FinalStateJets_{jet_radius:0.6f}_chg.dat"  

        n_wcut, temp_wcut = collect_jets_pTHatbin(fname_wcut, eta_cut, ptbins)
        n_ncut, temp_ncut = collect_jets_pTHatbin(fname_ncut, eta_cut, ptbins)
        #n_chgd, temp_chgd = collect_jets_pTHatbin(fname_chgd, eta_cut, ptbins)

        temp_ncut.multiply_by(temp_xsec)
        temp_wcut.multiply_by(temp_xsec)
        #temp_chgd.multiply_by(temp_xsec/n_chgd)
        nevt_ncut += n_ncut
        nevt_wcut += n_wcut
        jet_hist_ncut += temp_ncut
        jet_hist_wcut += temp_wcut
        #jet_hist_chgd += temp_chgd

    #scale = 1./nruns
    jet_hist_ncut.multiply_by(1./nevt_ncut)
    jet_hist_wcut.multiply_by(1./nevt_wcut)
    #jet_hist_chgd.multiply_by(scale)

    return jet_hist_ncut, jet_hist_wcut#, jet_hist_chgd 

def build_jet_spec(dirname, ptHatbins, nruns, eta_cut, ptbins, jet_radius):
    jhist_ncut = Histogram(ptbins) 
    jhist_wcut = Histogram(ptbins)   
    #jhist_chgd = Histogram(ptbins)
    num_pthatbins = len(ptHatbins)-1
    for ipT in range(num_pthatbins):
        ptmin, ptmax = ptHatbins[ipT], ptHatbins[ipT+1]
        #print(f"pT hat bin: {ptmin} --> {ptmax}")
        pthatdir = dirname + f"{ptmin}_{ptmax}_" + "{0:d}/"
        temp_ncut, temp_wcut = average_pTHat(pthatdir, nruns, eta_cut, ptbins, jet_radius)
        #, temp_chgd = average_pTHat(pthatdir, nruns, eta_cut, ptbins, jet_radius)

        jhist_ncut += temp_ncut
        jhist_wcut += temp_wcut
        #jhist_chgd += temp_chgd

    ## at this point all pT hat bins have been averaged and included
    ## can return.
    return jhist_ncut, jhist_wcut#, jhist_chgd

def save_to_file(fname, hist1, hist2):#, hist3):
   
    ptbins = hist1.get_xbins()

    #ncut,dncut = hist1.get_content()
    #wcut,dwcut = hist2.get_content()
    #chgd,dchgd = hist3.get_content()

    with open(fname, 'w') as f:
        f.write("ptmin,ptmax,Ncut,dNcut,Wcut,dWcut\n")#),Chg,dChg\n")
        for i in range(len(ptbins)-1):
            ptmin = ptbins[i]
            ptmax = ptbins[i+1]
            n, dn = hist1[i]#ncut[i], dncut[i]
            w, dw = hist2[i]#wcut[i], dwcut[i]
            #c, dc = chgd[i], dchgd[i]
            f.write(f"{ptmin:0.2f},{ptmax:0.2f},")
            f.write(f"{n:0.5e},{dn:0.5e},{w:0.5e},{dw:0.5e}\n")#,{c:0.5e},{dc:0.5e}\n")
    print(f"done writing {fname}")

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
    
    prop = obs.jet_spectra[index_of_system]
    pThats = obs.pT_hat_bins[index_of_system]
    num_runs = obs.nruns[tag]

    save_name = '/home/rmyazdi/scratch/JETSCAPE_local/jetscape_exec/gamma_proj/results/'
    if maxT != '':
        save_name += f'maxT_{maxT}/'
    if eloss != '':
        save_name += f'{eloss}/'
    save_name += f'{system}_'
    if centrality != '':
        save_name += f'{centrality}_'

    save_name += 'jet_spec_jet_rad_{radius}_{eta:0.2f}.csv'
    jet_radii = [0.2, 0.3, 0.4, 0.6, 0.8]
    print("JET SPECTRUM")
    for R in jet_radii:
        print(f"Start with Jet Radius {R}")
        ## first binning: just spectra for the large
        ## PseudoRapidity window.
        result = build_jet_spec(location, pThats, num_runs, prop['eta cut'], prop['spectrum bins'], R)
        ncutHist, wcutHist = result #, chgdHist = result

        ## Save these guys to file
        fname = save_name.format(radius=R, eta=prop['eta cut'])
        save_to_file(fname, ncutHist, wcutHist)#, chgdHist)

        del result, ncutHist, wcutHist, fname
        # Second Binning: eta_cut = 0.7 - R, so only apply for R < 0.7
        if R > 0.7:
            continue
        eta_wind = 0.7 - R
        fname = save_name.format(radius=R, eta=eta_wind)
        result = build_jet_spec(location, pThats, num_runs, eta_wind, prop['spectrum bins'], R)
        #ncutHist, wcutHist, chgdHist = result
        save_to_file(fname, *result)

