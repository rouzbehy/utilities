import os
import sys
from numpy import sqrt
import PDG
from Helpers import Histogram, Jet, Particle, get_xsection
import observables as obs
from typing import List, Dict, Tuple

def collect_fragmentation_funcs(fname: str, jet_radius: float, args: Dict) -> Tuple[int, Histogram, Histogram, Histogram]:
    """
    Collects data from a file and fills histograms with jet spectrum, fragmentation functions (FFs), and other relevant information.
    :param fname: The name of the file to read from.
    :param jet_radius: The radius of the jet.
    :param args: A dictionary containing various parameters.
    :return: The number of events, the jet spectrum histogram, and the FF histograms for z and pT.
    """
    jet_spec  = Histogram(args['spec bins'])
    ff_z      = Histogram(args['z bins'])
    ff_pT     = Histogram(args['pT bins'])
    nevt = 0
    ## cuts:
    pT_min_ff  = args['spec bins'][0]
    pT_max_ff  = args['spec bins'][-1]
    ff_eta_cut = args['eta cut']
    currjet = 0
    process_jet_ff = False
    for line in open(fname,'r'):
        line = line.split(',')
        len_line = len(line)
        if len_line == 1:
            nevt += 1
            continue
        elif len_line == 5:## it's a jet
            currjet = Jet(line)
            process_jet_ff = False
            if abs(currjet.eta) < ff_eta_cut and pT_min_ff <= currjet.pT < pT_max_ff:
                jet_spec.fill(currjet.pT)
                process_jet_ff = True
        elif len_line == 7: ## it's a jet constituent
            if not process_jet_ff:
                continue
            particle = Particle(line)
            if not (currjet.jetNum == particle.jetNum):
                print(f"Jet Number {currjet.jetNum} doesn't match the particle's jet number {particle.jetNum}")
                exit(-1)
            if particle.hadID not in PDG.PDG_Pool : continue
            if PDG.PDG_Charge[particle.hadID] == 0 : continue
            z = currjet.calc_z(particle)
            #if process_jet_ff:
            ff_z.fill(z)
            ff_pT.fill(particle.pT)
        else:
            continue
    return nevt, jet_spec, ff_z, ff_pT

def average_pTHat_bin_frag_func(dirname, num_runs, jet_radius, args):
    """
        Averages the fragmentation functions over a specific pT hat bin.
        :param dirname: The directory containing the files to read from.
        :param num_runs: The number of runs to average over.
        :param jet_radius: The radius of the jet.
        :param args: A dictionary containing various parameters.
        :return: The averaged histograms for the jet spectrum and FFs.
    """
    ncut_jet_spec  = Histogram(args['spec bins'])
    ncut_ff_z      = Histogram(args['z bins'])
    ncut_ff_pT     = Histogram(args['pT bins'])

    wcut_jet_spec  = Histogram(args['spec bins'])
    wcut_ff_z      = Histogram(args['z bins'])
    wcut_ff_pT     = Histogram(args['pT bins'])

    #chgd_jet_spec  = Histogram(args['spec bins'])
    #chgd_ff_z      = Histogram(args['z bins'])
    #chgd_ff_pT     = Histogram(args['pT bins'])

    nevt_ncut, nevt_wcut = 0, 0 #, nevt_chgd = 0,0,0
    for irun in range(num_runs):
        curr_dir = dirname.format(irun)
        temp_xsec = get_xsection(curr_dir)

        fname_wcut = curr_dir + f"FinalStateJets_{jet_radius:0.6f}_wcut.dat"
        fname_ncut = curr_dir + f"FinalStateJets_{jet_radius:0.6f}_ncut.dat"
        #fname_chgd = curr_dir + f"FinalStateJets_{jet_radius:0.6f}_chg.dat"

        n_ncut, tmp_ncut_jspec, tmp_ncut_ff_z, tmp_ncut_ff_pT = collect_fragmentation_funcs(fname_ncut, jet_radius, args)
        n_wcut, tmp_wcut_jspec, tmp_wcut_ff_z, tmp_wcut_ff_pT = collect_fragmentation_funcs(fname_wcut, jet_radius, args)
        #n_chgd, tmp_chgd_jspec, tmp_chgd_ff_z, tmp_chgd_ff_pT = collect_fragmentation_funcs(fname_chgd, jet_radius, args)

        tmp_scale_ncut = temp_xsec/n_ncut
        tmp_scale_wcut = temp_xsec/n_wcut
        #tmp_scale_chgd = temp_xsec/n_chgd

        tmp_ncut_jspec.multiply_by(tmp_scale_ncut)
        tmp_ncut_ff_z.multiply_by (tmp_scale_ncut)
        tmp_ncut_ff_pT.multiply_by(tmp_scale_ncut)

        tmp_wcut_jspec.multiply_by(tmp_scale_wcut)
        tmp_wcut_ff_z.multiply_by (tmp_scale_wcut)
        tmp_wcut_ff_pT.multiply_by(tmp_scale_wcut)

        #tmp_chgd_jspec.multiply_by(tmp_scale_chgd)
        #tmp_chgd_ff_z.multiply_by (tmp_scale_chgd)
        #tmp_chgd_ff_pT.multiply_by(tmp_scale_chgd)

        ncut_jet_spec += tmp_ncut_jspec
        ncut_ff_z     += tmp_ncut_ff_z
        ncut_ff_pT    += tmp_ncut_ff_pT

        wcut_jet_spec += tmp_wcut_jspec
        wcut_ff_z     += tmp_wcut_ff_z
        wcut_ff_pT    += tmp_wcut_ff_pT

        #chgd_jet_spec += tmp_chgd_jspec
        #chgd_ff_z     += tmp_chgd_ff_z
        #chgd_ff_pT    += tmp_chgd_ff_pT

    ## TODO use the number of events as weight
    scale = 1./num_runs
    ncut_jet_spec.multiply_by(scale)
    ncut_ff_z.multiply_by(scale)
    ncut_ff_pT.multiply_by(scale)

    wcut_jet_spec.multiply_by(scale)
    wcut_ff_z .multiply_by(scale)
    wcut_ff_pT.multiply_by(scale)

    #chgd_jet_spec.multiply_by(scale)
    #chgd_ff_z .multiply_by(scale)
    #chgd_ff_pT.multiply_by(scale)

    return ncut_jet_spec, ncut_ff_z, ncut_ff_pT,\
           wcut_jet_spec, wcut_ff_z, wcut_ff_pT
#chgd_jet_spec, chgd_ff_z, chgd_ff_pT

def build_fragmentation_functions(dirname, ptHatbins, nruns, jet_radius, args):

    ncut_jet_spec  = Histogram(args['spec bins'])
    ncut_ff_z      = Histogram(args['z bins'])
    ncut_ff_pT     = Histogram(args['pT bins'])

    wcut_jet_spec  = Histogram(args['spec bins'])
    wcut_ff_z      = Histogram(args['z bins'])
    wcut_ff_pT     = Histogram(args['pT bins'])

    #chgd_jet_spec  = Histogram(args['spec bins'])
    #chgd_ff_z      = Histogram(args['z bins'])
    #chgd_ff_pT     = Histogram(args['pT bins'])

    num_pthatbins = len(ptHatbins)-1
    for ipT in range(num_pthatbins):
        ptmin, ptmax = ptHatbins[ipT], ptHatbins[ipT+1]
        print(f"pT Hat bin: {ptmin} --> {ptmax}")
        pthatdir = dirname + f"{ptmin}_{ptmax}_" + "{0:d}/"
        hists = average_pTHat_bin_frag_func(pthatdir, nruns, jet_radius, args)

        ncut_jet_spec += hists[0]
        ncut_ff_z     += hists[1]
        ncut_ff_pT    += hists[2]

        wcut_jet_spec += hists[3]
        wcut_ff_z     += hists[4]
        wcut_ff_pT    += hists[5]

        #chgd_jet_spec += hists[6]
        #chgd_ff_z     += hists[7]
        #chgd_ff_pT    += hists[8]

    return ncut_jet_spec, ncut_ff_z, ncut_ff_pT,\
           wcut_jet_spec, wcut_ff_z, wcut_ff_pT
#chgd_jet_spec, chgd_ff_z, chgd_ff_pT

def save_to_file(fname, hist1, hist2, header):#, hist3, header):

    x = hist1.get_xbins()

    #ncut,dncut = hist1.get_content()
    #wcut,dwcut = hist2.get_content()
    #chgd,dchgd = hist3.get_content()

    with open(fname, 'w') as f:
        f.write(f"{header}\n")
        for i in range(len(x)-1):
            xmin = x[i]
            xmax = x[i+1]
            n, dn = hist1[i]#ncut[i], dncut[i]
            w, dw = hist2[i]#wcut[i], dwcut[i]
            #c, dc = chgd[i], dchgd[i]
            f.write(f"{xmin:0.4f},{xmax:0.4f},")
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
        maxT  = sys.argv[2]
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

    prop_LHC  = obs.jet_FF_LHC
    prop_RHIC = obs.jet_FF_RHIC
    pTHats = obs.pT_hat_bins[index_of_system]
    num_runs = obs.nruns[tag]

    save_name = '/home/rmyazdi/scratch/JETSCAPE_local/jetscape_exec/gamma_proj/results/'
    if maxT != '':
        save_name += f'maxT_{maxT}/'
    if eloss != '':
        save_name += f'{eloss}/'
    jet_radii = [0.2, 0.3, 0.4]#, 0.6, 0.8]

    for R in jet_radii:
        print(f"Start with Jet Radius {R}")
        ## first pass: LHC only
        if index_of_system != 0:
            result = build_fragmentation_functions(location, pTHats, num_runs, R, prop_LHC)
            ## save to file:
            header_jet_spec = "pTmin,pTmax,Ncut,dNcut,Wcut,dWcut"#,Chgd,dChgd"
            eta = prop_LHC['eta cut']
            fname1 = save_name + f'{system}_'
            if centrality != '':
                fname1 += f'{centrality}_'
            fname1 += f'LHC_jet_spec_FF_jet_rad_{R:0.1f}_{eta:0.2f}.csv'
            save_to_file(fname1, result[0], result[3], header_jet_spec)#, result[6], header_jet_spec)

            header_ff_z = "zmin,zmax,Ncut,dNcut,Wcut,dWcut"#,Chgd,dChgd"
            eta = prop_LHC['eta cut']
            fname2 = save_name + f'{system}_'
            if centrality != '':
                fname2 += f'{centrality}_'
            fname2 += f'LHC_FF_z_jet_rad_{R:0.1f}_{eta:0.2f}.csv'
            save_to_file(fname2, result[1], result[4],header_ff_z)#, result[7], header_ff_z)

            header_ff_pT = "pTmin,pTmax,Ncut,dNcut,Wcut,dWcut"#,Chgd,dChgd"
            eta = prop_LHC['eta cut']
            fname3 = save_name + f'{system}_'
            if centrality != '':
                fname3 += f'{centrality}_'
            fname3 += f'LHC_FF_pT_jet_rad_{R:0.1f}_{eta:0.2f}.csv'
            save_to_file(fname3, result[2], result[5], header_ff_pT)#, result[8], header_ff_pT)

        result = build_fragmentation_functions(location, pTHats, num_runs, R, prop_RHIC)
        ## save to file:
        header_jet_spec = "pTmin,pTmax,Ncut,dNcut,Wcut,dWcut"#,Chgd,dChgd"
        eta = prop_RHIC['eta cut']
        fname1 = save_name + f'{system}_'
        if centrality != '':
            fname1 += f'{centrality}_'
        fname1 += f'RHIC_jet_spec_FF_jet_rad_{R:0.1f}_{eta:0.2f}.csv'
        save_to_file(fname1, result[0], result[3],header_jet_spec)#, result[6], header_jet_spec)

        header_ff_z = "zmin,zmax,Ncut,dNcut,Wcut,dWcut"#,Chgd,dChgd"
        eta = prop_RHIC['eta cut']
        fname2 = save_name + f'{system}_'
        if centrality != '':
            fname2 += f'{centrality}_'
        fname2 += f'RHIC_FF_z_jet_rad_{R:0.1f}_{eta:0.2f}.csv'
        save_to_file(fname2, result[1], result[4], header_ff_z)#, result[7], header_ff_z)

        header_ff_pT = "pTmin,pTmax,Ncut,dNcut,Wcut,dWcut"#,Chgd,dChgd"
        eta = prop_RHIC['eta cut']
        fname3 = save_name + f'{system}_'
        if centrality != '':
            fname3 += f'{centrality}_'
        fname3 += f'RHIC_FF_pT_jet_rad_{R:0.1f}_{eta:0.2f}.csv'
        save_to_file(fname3, result[2], result[5], header_ff_pT)#, result[8], header_ff_pT)

    print("Done with the jet fragmentation functions")
