import os
import sys
from numpy import sqrt, cos, sin
import PDG
from Helpers import Histogram, Jet, Particle,get_xsection
import observables as obs

def collect_jet_shapes(fname, jetR, rbins, eta_min, eta_max, pTmin, pTmax):
    jet_shape = Histogram(rbins)
    nevt = 0
    currjet = 0
    process_jet_shape = False
    num_jets = 0
    for line in open(fname,'r'):
        line = line.split(',')
        len_line = len(line)
        if len_line == 1:
            nevt += 1
            continue
        elif len_line == 5:## it's a jet
            currjet = Jet(line)
            process_jet_shape = False
            if eta_min < abs(currjet.eta) < eta_max and pTmin < currjet.pT < pTmax:
                process_jet_shape = True
                num_jets += 1
        elif len_line == 7: ## it's a jet constituent
            if not process_jet_shape:
                continue
            particle = Particle(line)
            ## make sure the particle here corresponds to the same jet
            if not (currjet.jetNum == particle.jetNum):
                print(f"Jet {currjet.jetNum} doesn't match particle's jnumber {particle.jetNum}")
                exit(-1)
            if particle.hadID not in PDG.PDG_Pool : continue
            if PDG.PDG_Charge[particle.hadID] == 0 : continue
            ## calculate the r
            r = currjet.calc_r(particle)
            jet_shape.fill(r, particle.pT/currjet.pT)
        else:
            continue
    return nevt, num_jets, jet_shape

def average_pTHat(dirname, nruns, jetR, rbins, eta_min, eta_max, pTmin, pTmax):

    jet_hist_wcut = Histogram(rbins)
    jet_hist_ncut = Histogram(rbins)
    ##jet_hist_chgd = Histogram(rbins)
    nevents_ncut, nevents_wcut = 0, 0
    njets_ncut, njets_wcut = 0, 0
    for irun in range(nruns):
        curr_dir = dirname.format(irun)
        temp_xsec = get_xsection(curr_dir)

        fname_wcut = curr_dir + f"FinalStateJets_{jetR:0.6f}_wcut.dat" 
        fname_ncut = curr_dir + f"FinalStateJets_{jetR:0.6f}_ncut.dat" 
        ##fname_chgd = curr_dir + f"FinalStateJets_{jetR:0.6f}_chg.dat"  

        n_wcut, nj_wcut, temp_wcut = collect_jet_shapes(fname_wcut, jetR, rbins, eta_min, eta_max, pTmin, pTmax)
        n_ncut, nj_ncut, temp_ncut = collect_jet_shapes(fname_ncut, jetR, rbins, eta_min, eta_max, pTmin, pTmax)
        ##n_chgd, nj_chgd, temp_chgd = collect_jet_shapes(fname_chgd, jetR, rbins, eta_min, eta_max, pTmin, pTmax)
        ## TODO : number of jets should not be divided here
        nevents_ncut += n_ncut
        nevents_wcut += n_wcut
        njets_ncut += nj_ncut
        njets_wcut += nj_wcut

        temp_ncut.multiply_by(temp_xsec)
        temp_wcut.multiply_by(temp_xsec)
        ##temp_chgd.multiply_by(temp_xsec/(n_chgd*nj_chgd))

        jet_hist_ncut += temp_ncut
        jet_hist_wcut += temp_wcut
        ##jet_hist_chgd += temp_chgd
    ## TODO divide by number of events
    #jet_hist_ncut.multiply_by(1./nruns)
    #jet_hist_wcut.multiply_by(1./nruns)
    ##jet_hist_chgd.multiply_by(1./nruns)
    return nevents_ncut, njets_ncut ,jet_hist_ncut,\
           nevents_wcut, njets_wcut, jet_hist_wcut#,\
           #, jet_hist_chgd

def build_jet_shape(dirname, pTHatbins, nruns, jetR, rbins, eta_min, eta_max, pTmin, pTmax, debug=False):
    
    jhist_ncut = Histogram(rbins) 
    jhist_wcut = Histogram(rbins)   
    ##jhist_chgd = Histogram(rbins)
    num_jets_ncut, num_jets_wcut = 0, 0
    num_evts_ncut, num_evts_wcut = 0, 0
    debug_dirname = "./debug/"
    if not os.path.exists(debug_dirname):
        os.system(f"mkdir -p {debug_dirname}")
    num_pthatbins = len(pTHatbins)-1
    for ipT in range(num_pthatbins):
        ptHatmin, ptHatmax = pTHatbins[ipT], pTHatbins[ipT+1]
        #print(f"Dealing with pTHat bin: {ptHatmin} --> {ptHatmax}")
        pthatdir = dirname + f"{ptHatmin}_{ptHatmax}_" + "{0:d}/"
        res = average_pTHat(pthatdir, nruns, jetR, rbins, eta_min, eta_max, pTmin, pTmax)
        (nevts_ncut, njets_ncut, temp_ncut, nevts_wcut, njets_wcut, temp_wcut) = res
        #, temp_chgd = average_pTHat(pthatdir, nruns, jetR, rbins, eta_min, eta_max, pTmin, pTmax)
        num_evts_ncut += nevts_ncut 
        num_jets_ncut += njets_ncut 
        num_evts_wcut += nevts_wcut 
        num_jets_wcut += njets_wcut 
        if debug:
            fname = debug_dirname + f"jet_shape_pTHat_{ipT}_R_{jetR}_{eta_min:0.1f}_{eta_max:0.1f}.csv"
            save_to_file(fname, pTmin, pTmax,temp_ncut, temp_wcut, nevts_ncut, njets_ncut, nevts_wcut, njets_wcut)
        jhist_ncut += temp_ncut
        jhist_wcut += temp_wcut
        #jhist_chgd += temp_chgd

    return jhist_ncut, num_evts_ncut, num_jets_ncut,\
           jhist_wcut, num_evts_wcut, num_jets_wcut
           #, jhist_chgd


def save_to_file(fname, ptmin, ptmax, hist1, hist2, nevt1, njets1, nevt2, njets2):#, hist3):
    xbins = hist1.get_xbins()

    with open(fname, 'w') as f:
        f.write(f"#ptmin {ptmin} ptmax {ptmax} ")
        f.write(f"nevt_ncut {nevt1} njet_ncut {njets1} ")
        f.write(f"nevt_wcut {nevt2} njet_wcut {njets2}\n")
        f.write("rmin,rmax,ncut,dncut,wcut,dwcut\n")#,chgd,dchgd\n")
        #f.write("ptmin,ptmax,Ncut,dNcut,Wcut,dWcut,Chg,dChg\n")
        for i in range(len(xbins)-1):
            xmin = xbins[i]
            xmax = xbins[i+1]
            n, dn = hist1[i]#ncut[i], dncut[i]
            w, dw = hist2[i]#wcut[i], dwcut[i]
            #c, dc = chgd[i], dchgd[i]
            f.write(f"{xmin:0.2f},{xmax:0.2f},")
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
        eloss       = sys.argv[2]
        centrality  = sys.argv[3]
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
    
    prop_LHC = obs.jet_shape_LHC
    prop_RHIC = obs.jet_shape_RHIC
    pTHats = obs.pT_hat_bins[index_of_system]
    num_runs = obs.nruns[tag]

    save_name = '/home/rmyazdi/scratch/JETSCAPE_local/jetscape_exec/gamma_proj/results/'
    if maxT != '':
        save_name += f'maxT_{maxT}/'

    if eloss != '':
        save_name += f'{eloss}/'
    save_name += f'{system}_'
    if centrality != '':
        save_name += f'{centrality}_'
    save_name += 'jet_shape_{cutType}_cuts_jet_rad_{radius}_{etamin:0.2f}_{etamax:0.2f}.csv'
    jet_radii = [0.2, 0.3, 0.4]#, 0.6, 0.8]
    print("PROCESS JET SHAPE")
    for R in jet_radii:
        print(f"Start with Jet Radius {R}")
        ## first binning: just spectra for the large
        ## PseudoRapidity window.
        result = build_jet_shape(location, pTHats, num_runs, R, 
                                 prop_LHC['r bins'], prop_LHC['eta min'], 
                                 prop_LHC['eta max'], prop_LHC['pT min'], 
                                 prop_LHC['pT max'], debug=True)
        jhist_ncut, num_evts_ncut, num_jets_ncut, jhist_wcut, num_evts_wcut, num_jets_wcut = result
        #ncutHist, wcutHist, chgdHist = result

        ## Save these guys to file
        fname = save_name.format(radius=R, 
                                etamin=prop_LHC['eta min'], etamax=prop_LHC['eta max'],
                                cutType='LHC')

        save_to_file(fname, prop_LHC['pT min'], prop_LHC['pT max'], \
                    jhist_ncut, jhist_wcut, num_evts_ncut, num_jets_ncut, \
                    num_evts_wcut, num_jets_wcut)#, chgdHist)

        #del result, ncutHist, wcutHist#, chgdHist 

        #fname = save_name.format(radius=R, 
        #                        etamin=prop_RHIC['eta min'], etamax=prop_RHIC['eta max'],
        #                        cutType='RHIC')

        #result = build_jet_shape(location, pTHats, num_runs, R, 
        #                         prop_RHIC['r bins'] , prop_RHIC['eta min'], 
        #                         prop_RHIC['eta max'], prop_RHIC['pT min'], 
        #                         prop_RHIC['pT max'])

        #ncutHist, wcutHist = result #, chgdHist = result
        #save_to_file(fname, prop_RHIC['pT min'], prop_RHIC['pT max'], ncutHist, wcutHist)#, chgdHist)
