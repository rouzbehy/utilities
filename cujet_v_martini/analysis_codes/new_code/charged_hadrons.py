import sys
import os

# from numpy import array, sqrt
import PDG
from Helpers import Histogram, get_xsection, Particle
import observables as obs
from typing import Type


def get_hadrons_single_file(fname: str, hist: Histogram, ch_had_eta_cut: float):
    # infile     = directory + "FinalStateHadrons.txt"
    # if not os.path.exists(infile):
    #    print(f"Cannot open the input file {infile}")
    #    return False
    nevt = 0
    for line in open(fname, "r"):
        if "#Event_ID" in line:
            nevt += 1
            continue
        line = line.split(",")
        if not len(line) == 7:
            continue
        (index, pstat, pid, en, pT, phi, eta) = line
        pstat = int(pstat)
        weight = 1 if pstat >= 0 else -1
        pid = abs(int(pid))
        en = float(en)
        pT = float(pT)
        phi = float(phi)
        eta = float(eta)
        if pid not in PDG.PDG_Pool:
            continue
        if PDG.PDG_Charge[pid] == 0:
            continue
        if abs(eta) < ch_had_eta_cut:
            hist.fill(pT, weight)
    return nevt, hist


def collect_hadrons_pTHatbin(data_loc, pTmin, pTmax, num_runs, pTbins, ch_had_eta_cut):
    """
    collect_hadrons_pTHatbin(data_loc, pTmin, pTmax, num_runs, hist):
        - data_loc: where to look
        - fname: file to process
        - hist: the histogram to fill
    return:
        histogram
    """
    n_ch = []
    print(f"working on {pTmin}, {pTmax}")
    nevts = 0
    for irun in range(num_runs):
        directory = data_loc + f"{pTmin:d}_{pTmax:d}_{irun}/"
        fname = directory + "FinalStateHadrons.txt"
        temp_hist = Histogram(pTbins)
        temp_nevt, temp_hist = get_hadrons_single_file(fname, temp_hist, ch_had_eta_cut)
        cross_section = get_xsection(directory)
        ## scale the histogram with the cross section
        temp_hist.multiply_by(cross_section)
        nevts += temp_nevt
        n_ch.append(temp_hist)

    spec = Histogram(pTbins)
    for hist in n_ch:
        spec += hist
    spec.multiply_by(1.0 / nevts)
    return spec


def build_charged_spec(data_loc, pTHatbins, num_runs, pTbins, ch_had_eta_cut):
    result_hist = Histogram(pTbins)
    num_pTHat_bins = len(pTHatbins) - 1
    for ibin in range(num_pTHat_bins):
        pTmin, pTmax = pTHatbins[ibin], pTHatbins[ibin + 1]
        tmp_hist = collect_hadrons_pTHatbin(
            data_loc, pTmin, pTmax, num_runs, pTbins, ch_had_eta_cut
        )
        result_hist = result_hist + tmp_hist

    return result_hist


def save_charged_spec(fname, hist):
    # spectrum, errs = hist.get_content()
    x = hist.get_xbins()
    with open(fname, "w") as out:
        out.write("pTmin,pTmax,N,dN\n")
        num_bins = hist.nbins
        for ibin in range(num_bins):
            x1, x2 = x[ibin], x[ibin + 1]
            y, dy = hist[ibin]  # spectrum[ibin], errs[ibin]
            out.write(f"{x1:0.2f},{x2:0.2f},{y:0.5e},{dy:0.5e}\n")
    print(f"Finished writing {fname}")


if __name__ == "__main__":
    ## get the system
    system = sys.argv[1]
    eloss = ""
    centrality = ""
    maxT = ""

    location = "/home/rmyazdi/scratch/JETSCAPE_local/jetscape_exec/gamma_proj/"
    if "pp" not in system:
        eloss = sys.argv[2]
        centrality = sys.argv[3]
        location += f"HIC/{eloss}/{system}_{centrality}/"
    else:
        maxT = sys.argv[2]
        location += f"pp/results/maxT_{maxT}/{system}/"
    index_of_system = -1
    if "200" in system:
        index_of_system = 0
    elif "2760" in system:
        index_of_system = 1
    else:
        index_of_system = 2

    tag = "new"
    if system == "PbPb2760" and centrality in ["00-05", "20-30", "30-40"]:
        if eloss != "cujet":
            tag = "old"
    prop = obs.charged_hadrons[index_of_system]
    pThats = obs.pT_hat_bins[index_of_system]
    num_runs = obs.nruns[tag]

    hist = build_charged_spec(
        location, pThats, num_runs, prop["spectrum bins"], prop["eta cut"]
    )

    save_name = "/home/rmyazdi/scratch/JETSCAPE_local/jetscape_exec/gamma_proj/results/"
    if maxT != "":
        save_name += f"maxT_{maxT}/"
    if eloss != "":
        save_name += f"{eloss}/"

    if not os.path.exists(save_name):
        os.system(f"mkdir -p {save_name}")
    eta = prop["eta cut"]
    save_name += f"{system}_"
    if centrality != "":
        save_name += f"{centrality}_"
    save_name += f"charged_hadrons_eta_cut_{eta}.csv"
    save_charged_spec(save_name, hist)
