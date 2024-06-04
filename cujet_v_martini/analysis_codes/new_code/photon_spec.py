import sys
import os
from Helpers import Histogram, get_xsection, Particle
import observables as obs

def collect_photons_pTHatbin(fname, photon_eta_cut, ptbins):
    """
        def collect_photons_pTHatbin(fname, photon_eta_cut, ptbins):
            read the file with _fname_ and populate 
            the `photon_hist` object with photons within
            the given `eta_cut` window
    """
    nevt = 0
    pt = 0.
    stat = 0
    photon_hist_conv = Histogram(ptbins)
    photon_hist_brem = Histogram(ptbins)
    photon_hist_othr = Histogram(ptbins)
    for line in open(fname,'r'):
        line = line.split(',')
        len_line = len(line)
        if "#Event" in line[0]:
            nevt += 1
            continue
        elif len_line == 8:## it's a particle
            #"P\ti,ipart,p.pid(),p.pstat(),p.e(),p.pt(),p.phi(),p.eta()"
            if int(line[2]) == 22: ## a photon:
                stat = int(line[3])
                pt = float(line[5])
                eta = float(line[7])
                if abs(eta) > photon_eta_cut:
                    continue
                if stat == 99:
                    photon_hist_conv.fill(pt)
                elif stat == 100:
                    photon_hist_brem.fill(pt)
                elif stat != 91:## don't want decay photons from pi^0 or eta
                    ##print(f"Have a non jet-medium photon with pT: {pt}")
                    photon_hist_othr.fill(pt)
                else:
                    continue
        else:
            continue
    return nevt, photon_hist_conv, photon_hist_brem, photon_hist_othr

def average_pTHat(dirname, nruns, eta_cut, ptbins):
    photon_hist_conv = Histogram(ptbins)
    photon_hist_brem = Histogram(ptbins)
    photon_hist_othr = Histogram(ptbins)
    nevt = 0
    for irun in range(nruns):
        curr_dir = dirname.format(irun)
        temp_xsec = get_xsection(curr_dir)
        fname = curr_dir + f"FinalStateParticles.txt"
        result = collect_photons_pTHatbin(fname, eta_cut, ptbins)
        (tmp_nevt, tmp_hist_conv, tmp_hist_brem, tmp_hist_othr) = result
        nevt += tmp_nevt
        tmp_hist_conv.multiply_by(temp_xsec)
        tmp_hist_brem.multiply_by(temp_xsec)
        tmp_hist_othr.multiply_by(temp_xsec)
        photon_hist_conv += tmp_hist_conv
        photon_hist_brem += tmp_hist_brem
        photon_hist_othr += tmp_hist_othr

    scale=1./nevt
    photon_hist_conv.multiply_by(scale)
    photon_hist_brem.multiply_by(scale)
    photon_hist_othr.multiply_by(scale)

    return photon_hist_conv, photon_hist_brem, photon_hist_othr

def build_photon_spec(dirname, ptHatbins, nruns, eta_cut, ptbins):
    photon_hist_conv = Histogram(ptbins)
    photon_hist_brem = Histogram(ptbins)
    photon_hist_othr = Histogram(ptbins)
    num_pthatbins = len(ptHatbins)-1
    for ipT in range(num_pthatbins):
        ptmin, ptmax = ptHatbins[ipT], ptHatbins[ipT+1]
        pthatdir = dirname + f"{ptmin}_{ptmax}_" + "{0:d}/"
        temp_conv, temp_brem, temp_othr = average_pTHat(pthatdir, nruns, eta_cut, ptbins)
        photon_hist_conv += temp_conv
        photon_hist_brem += temp_brem
        photon_hist_othr += temp_othr

    return photon_hist_conv,photon_hist_brem,photon_hist_othr

def save_to_file(fname, hist1, hist2, hist3):
   
    ptbins = hist1.get_xbins()

    #conv, dconv = hist1.get_content()
    #brem, dbrem = hist2.get_content()
    #othr, dothr = hist3.get_content()

    with open(fname, 'w') as f:
        f.write("ptmin,ptmax,conv,dconv,brem,dbrem,othr,dothr\n")
        for i in range(len(ptbins)-1):
            ptmin = ptbins[i]
            ptmax = ptbins[i+1]
            v1, dv1 = hist1[i]#conv
            v2, dv2 = hist2[i]#brem
            v3, dv3 = hist3[i]#othr
            f.write(f"{ptmin:0.2f},{ptmax:0.2f},")
            f.write(f"{v1:0.5e},{dv1:0.5e},{v2:0.5e},{dv2:0.5e},")
            f.write(f"{v3:0.5e},{dv3:0.5e}\n")
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
        if eloss == 'martini':
            tag = 'old'
    
    prop = obs.photon_specs[index_of_system]
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

    save_name += 'photon_spec_{eta:0.2f}.csv'

    result = build_photon_spec(location, pThats, num_runs, prop['eta cut'], prop['spectrum bins'])
    Hist_Conv, Hist_Brem, Hist_Othr = result

    ## Save these guys to file
    fname = save_name.format(eta=prop['eta cut'])
    save_to_file(fname, Hist_Conv, Hist_Brem, Hist_Othr)
