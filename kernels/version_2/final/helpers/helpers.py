import pandas as pd
import numpy as np
import helpers.dictionaries as my_dicts
import helpers.jetDicts as my_jdicts
from typing import List


def get_experimental_charged_RAA(pT_low_cut: float, pT_high_cut: float):
    exp_loc = "/Users/rmyazdi/Documents/research/"
    exp_loc += "KERNELS_NLO_NP_PART2/exp_data/charged/"
    cms_cents = {1: "0_5", 2: "5_10", 3: "10_30"}
    cms = {
        c: pd.read_csv(exp_loc + f"CMS/Table5_{i}.csv", comment="#").rename(
            columns=my_dicts.colnames_raa_cms
        )
        for i, c in cms_cents.items()
    }

    alice_cents = {"00_05": "0_5", "05_10": "5_10", "10_20": "10_20"}
    alice = {
        c: pd.read_csv(
            exp_loc + f"ALICE/PbPb_RAA_cent_{i}_eta_0p8.csv", comment="#"
        ).rename(columns=my_dicts.colnames_raa_alice)
        for i, c in alice_cents.items()
    }

    atlas_1 = pd.read_csv(exp_loc + "ATLAS/Table49.csv", comment="#").rename(
        columns=my_dicts.colnames_raa_atlas
    )
    atlas_2_cents = {33: "0_5", 34: "5_10", 35: "10_20"}
    atlas_2 = {
        c: pd.read_csv(exp_loc + f"ATLAS/Table{i}.csv", comment="#").rename(
            columns=my_dicts.colnames_raa_atlas
        )
        for i, c in atlas_2_cents.items()
    }

    cms = {c: itm[itm["x"].between(pT_low_cut, pT_high_cut)] for c, itm in cms.items()}
    alice = {
        c: itm[itm["x"].between(pT_low_cut, pT_high_cut)] for c, itm in alice.items()
    }

    atlas_1 = atlas_1[atlas_1["x"].between(pT_low_cut, pT_high_cut)]
    atlas_2 = {
        c: itm[itm["x"].between(pT_low_cut, pT_high_cut)] for c, itm in atlas_2.items()
    }

    return {"cms": cms, "alice": alice, "atlas_1": atlas_1, "atlas_2": atlas_2}


def get_charged_RAA(prefix: str, run_type: str, pT_low_cut: float, pT_high_cut: float):

    ## read in the calculated pp hadron spec for the new formation calculations MARTINI
    pp_loc, aa_loc = prefix, prefix
    if run_type == "with formation time":
        pp_loc += "pp_2p76/chgd_spec.csv"
        aa_loc += "final_runs/rset_{rset}/2p76/{cent}/hadron_spectra.csv"
    else:
        pp_loc += "v2/proton_proton/data/hadron_spectra.csv"
        aa_loc += "v2/production_v2/martini_results/"
        aa_loc += "final_PbPb_2p76/rset_{rset}/cent_{cent}/hadron_spectra.csv"

    pp_data = pd.read_csv(pp_loc, comment="#")
    pp_data["pT"] = 0.5 * (pp_data["pTmin"] + pp_data["pTmax"])
    pp_data["dpT"] = pp_data["pTmax"] - pp_data["pTmin"]
    pp_data = pp_data[pp_data["pT"].between(pT_low_cut, pT_high_cut)]

    ## Read in the AA charged hadron results. First the new calculations with formation time
    martini = {}
    for r in [1, 2, 3]:
        martini[r] = {}
        for c in ["0_5", "5_10", "10_20"]:
            dat = pd.read_csv(aa_loc.format(cent=c, rset=r), comment="#")
            dat["pT"] = 0.5 * (dat["pTmin"] + dat["pTmax"])
            dat["dpT"] = dat["pTmax"] - dat["pTmin"]
            martini[r][c] = dat[dat["pT"].between(pT_low_cut, pT_high_cut)]

    xp, yp, dyp = pp_data["pT"], pp_data["N"], pp_data["dN"]
    RAA = {}
    for rset in martini:
        RAA[rset] = {}
        for cent in martini[rset]:
            AA = martini[rset][cent]
            ya, dya = AA["N"], AA["dN"]
            raa = ya / yp
            draa = raa * np.sqrt((dyp / yp) ** 2 + (dya / ya) ** 2)
            RAA[rset][cent] = pd.DataFrame({"pT": xp, "RAA": raa, "dRAA": draa})

    return RAA


def get_experimental_jet_RAA(pT_low_cut: float, pT_high_cut: float):
    ## read in the experimental RAA:
    prefix = "/Users/rmyazdi/Documents/research/"
    prefix += "KERNELS_NLO_NP_PART2/exp_data/jets/"
    cms = {}
    cents = {"00-05": "0_5", "05-10": "5_10", "10-30": "10_30"}
    radii = ["0p2", "0p3", "0p4"]
    for cent, label in cents.items():
        cms[label] = {}
        for R in radii:
            cms[label][R] = pd.read_csv(
                prefix + f"CMS/jetscape/cent_{cent}_R_{R}.csv",
                comment="#",
            ).rename(columns=my_jdicts.colnames_CMS_RAA)
    return cms


def get_jet_RAA(prefix: str, run_type: str, pT_low_cut: float, pT_high_cut: float):
    #### read in the calculated pp charged hadron spec
    pp_loc, aa_loc = prefix, prefix
    if run_type == "with formation time":
        pp_loc += "pp_2p76/jet_spec.csv"
        aa_loc += "final_runs/rset_{r}/2p76/{c}/jet_spectra.csv"
    else:
        pp_loc += "v2/proton_proton/data/jet_spectra.csv"
        aa_loc += "v2/production_v2/martini_results/final_PbPb_2p76"
        aa_loc += "/rset_{r}/cent_{c}/jet_spectra.csv"
    ## read in the proton-proton jet spectrum
    pp_data = pd.read_csv(pp_loc, comment="#")
    pp_data["pT"] = 0.5 * (pp_data["pTmin"] + pp_data["pTmax"])
    pp_data["dpT"] = pp_data["pTmax"] - pp_data["pTmin"]
    pp_data = pp_data[pp_data["pT"].between(pT_low_cut, pT_high_cut)]

    ## read in the Pb-Pb jet spectrum for each rate set, centrality and
    ## cone radius
    jet_RAA = {}
    for r in [1, 2, 3]:
        jet_RAA[r] = {}
        for c in ["0_5", "5_10", "10_20"]:
            jet_RAA[r][c] = {}
            raw = pd.read_csv(aa_loc.format(r=r, c=c), comment="#")
            raw["pT"] = 0.5 * (raw["pTmin"] + raw["pTmax"])
            raw = raw[raw["pT"].between(pT_low_cut, pT_high_cut)]
            for radius in ["0p2", "0p3", "0p4"]:
                xp, yp, dyp = (
                    pp_data["pT"],
                    pp_data["N" + radius],
                    pp_data["dN" + radius],
                )
                xa, ya, dya = raw["pT"], raw["N" + radius], raw["dN" + radius]
                raa = ya / yp
                draa = raa * np.sqrt((dyp / yp) ** 2 + (dya / ya) ** 2)
                jet_RAA[r][c][radius] = pd.DataFrame(
                    {"pT": xp, "raa": raa, "draa": draa}
                )

    return jet_RAA


def get_experimental_shape_ratio():
    ## read in the experimental results
    data_fname = "/Users/rmyazdi/Documents/research/"
    data_fname += "KERNELS_NLO_NP_PART2/exp_data/jet_shape/"
    data_fname += "CMS_PbPb_2760_Jet_ShapeRatio_{c}_pT=100-INF_R=0p3_eta=0p3-2p0.dat"
    cms_data = {}
    cents = {"00-10": "0_10", "10-30": "10_30"}
    for c in ["00-10", "10-30"]:
        tmp = np.loadtxt(
            data_fname.format(c=c), comments="#", unpack=True, delimiter="\t"
        )
        cms_data[cents[c]] = pd.DataFrame(
            {"x": tmp[0], "y": tmp[1], "dx": tmp[2], "dy": tmp[3]}
        )
    return cms_data


def _get_aa_jet_shape(fname: str):
    tmp = pd.read_csv(fname, comment="#")
    njet = 1
    with open(fname, "r") as f:
        line = f.readline()
        line = line.split(" ")[-1]
        njet = float(line)
    dat = tmp[tmp["rmax"] < 0.31]
    delta_r = tmp["rmax"] - tmp["rmin"]
    r = 0.5 * (tmp["rmax"] + tmp["rmin"])
    rho = tmp["N"] / (njet)
    drho = tmp["dN"] / (njet)
    norm = sum(rho.to_list())
    rho_normed = rho / (delta_r * norm)
    drho_normed = drho / (delta_r * norm)
    return pd.DataFrame({"rho": rho_normed, "drho": drho_normed, "r": r, "dr": delta_r})


def get_jet_shape_ratios(prefix: str, run_type: str):

    aa_loc, pp_loc = prefix, prefix
    if run_type == "with formation time":
        aa_loc += "final_runs/rset_{r}/2p76/{c}/jet_shape.csv"
        pp_loc += "pp_2p76/jet_shape.csv"
    else:
        aa_loc += "v2/production_v2/martini_results/final_PbPb_2p76"
        aa_loc += "/rset_{r}/cent_{c}/jet_shape.csv"
        pp_loc += "v2/proton_proton/data/jet_shape.csv"

    # get proton-proton baseline
    calcs = pd.read_csv(pp_loc, comment="#")
    njet = 1
    with open(pp_loc, "r") as f:
        line = f.readline()
        line = line.split(" ")[-1]
        njet = float(line)

    dat = calcs[calcs["rmax"] < 0.31]
    delta_r = dat["rmax"] - dat["rmin"]
    r = 0.5 * (dat["rmax"] + dat["rmin"])
    rho = dat["N"] / (njet)
    drho = dat["dN"] / (njet)
    norm = sum(rho.to_list())  #
    rho_normed = rho / (delta_r * norm)
    drho_normed = drho / (delta_r * norm)

    pp_data = pd.DataFrame(
        {"r": r, "dr": delta_r / 2, "rho": rho_normed, "drho": drho_normed}
    )

    ## Compute the jet shape ratio for each rate set and the two
    ## centralities of interest:
    result = {}
    for r in [1, 2, 3]:
        result[r] = {}
        for c in ["0_5", "5_10", "10_20"]:
            aa = _get_aa_jet_shape(aa_loc.format(r=r, c=c))
            ratio = aa["rho"] / pp_data["rho"]
            dratio = ratio * np.sqrt(
                (aa["drho"] / aa["rho"]) ** 2 + (pp_data["drho"] / pp_data["rho"]) ** 2
            )
            result[r][c] = pd.DataFrame(
                {
                    "ratio": ratio,
                    "dratio": dratio,
                    "r": pp_data["r"],
                    "dr": pp_data["dr"],
                }
            )

        ## create the 0-10 centrality class by averaging
        ## the 0-5 and 5-10
        rho_normed = 0.5 * (result[r]["0_5"]["ratio"] + result[r]["5_10"]["ratio"])
        # scale = sum(rho_normed)
        drho_normed = 0.5 * np.sqrt(
            result[r]["0_5"]["dratio"] ** 2 + result[r]["5_10"]["dratio"] ** 2
        )
        result[r]["0_10"] = pd.DataFrame(
            {
                "ratio": rho_normed,
                "dratio": drho_normed,
                "r": result[r]["0_5"]["r"],
                "dr": result[r]["0_5"]["dr"],
            }
        )

    return result


def get_experimenta_FF_ratio():
    exp_loc = "/Users/rmyazdi/Documents/research/"
    exp_loc += "KERNELS_NLO_NP_PART2/exp_data/jet_fragmentation/ATLAS/"
    ## read in ATLAS results
    data_pT = pd.read_csv(
        exp_loc + "HEPData-ins1511869-v1-csv/Table9.csv", comment="#"
    ).rename(columns=my_jdicts.colnames_ATLAS_FF_DPT)
    data_z = pd.read_csv(
        exp_loc + "HEPData-ins1511869-v1-csv/Table25.csv", comment="#"
    ).rename(columns=my_jdicts.colnames_ATLAS_FF_DZ)
    return data_pT, data_z


def get_pp_jet_ff(
    fname_pt: str, fname_z: str, pTmin: float, pTmax: float, zmin: float, zmax: float
):
    tmp_DPT = pd.read_csv(fname_pt, comment="#")
    tmp_DZ = pd.read_csv(fname_z, comment="#")
    tmp_DPT = tmp_DPT[tmp_DPT["pTmax"] < pTmax]
    tmp_DPT = tmp_DPT[tmp_DPT["pTmin"] > pTmin]
    tmp_DZ = tmp_DZ[tmp_DZ["zmax"] < zmax]
    tmp_DZ = tmp_DZ[tmp_DZ["zmin"] > zmin]
    njets_pp, nevts_pp = 0, 0
    with open(fname_pt, "r") as f:
        line = f.readline()
        line = line.split(" ")
        njets_pp = float(line[1])
        nevts_pp = float(line[-1])

    dpT = tmp_DPT["pTmax"] - tmp_DPT["pTmin"]
    pT = 0.5 * (tmp_DPT["pTmax"] + tmp_DPT["pTmin"])
    scale = njets_pp
    tmpy = tmp_DPT["N"] / (dpT * scale)
    dtmpy = tmp_DPT["dN"] / (dpT * scale)

    dz = tmp_DZ["zmax"] - tmp_DZ["zmin"]
    z = 0.5 * (tmp_DZ["zmax"] + tmp_DZ["zmin"])
    tmpyz = tmp_DZ["N"] / (dz * scale)
    dtmpyz = tmp_DZ["dN"] / (dz * scale)

    pp_pT = pd.DataFrame({"pT": pT, "dpT": dpT, "FF(pT)": tmpy, "dFF(pT)": dtmpy})
    pp_z = pd.DataFrame({"z": z, "dz": dz, "FF(z)": tmpyz, "dFF(z)": dtmpyz})
    return pp_pT, pp_z


def get_aa_jet_ff(
    fname_pt: str, fname_z: str, pTmin: float, pTmax: float, zmin: float, zmax: float
):
    tmp_DPT = pd.read_csv(fname_pt, comment="#")
    tmp_DZ = pd.read_csv(fname_z, comment="#")
    tmp_DPT = tmp_DPT[tmp_DPT["pTmax"] < pTmax]
    tmp_DPT = tmp_DPT[tmp_DPT["pTmin"] > pTmin]
    tmp_DZ = tmp_DZ[tmp_DZ["zmax"] < zmax]
    tmp_DZ = tmp_DZ[tmp_DZ["zmin"] > zmin]
    njets_pp, nevts_pp = 0, 0
    with open(fname_pt, "r") as f:
        line = f.readline()
        line = line.split(" ")
        njets_pp = float(line[1])
        nevts_pp = float(line[-1])

    dpT = tmp_DPT["pTmax"] - tmp_DPT["pTmin"]
    pT = 0.5 * (tmp_DPT["pTmax"] + tmp_DPT["pTmin"])
    scale = njets_pp
    tmpy = tmp_DPT["N"] / (dpT * scale)
    dtmpy = tmp_DPT["dN"] / (dpT * scale)

    dz = tmp_DZ["zmax"] - tmp_DZ["zmin"]
    z = 0.5 * (tmp_DZ["zmax"] + tmp_DZ["zmin"])
    tmpyz = tmp_DZ["N"] / (dz * scale)
    dtmpyz = tmp_DZ["dN"] / (dz * scale)
    result_pT = pd.DataFrame({"pT": pT, "dpT": dpT, "FF(pT)": tmpy, "dFF(pT)": dtmpy})
    result_z = pd.DataFrame({"z": z, "dz": dz, "FF(z)": tmpyz, "dFF(z)": dtmpyz})

    return result_pT, result_z


def get_jet_FF_ratios(
    prefix: str, run_type: str, pTmin: float, pTmax: float, zmin: float, zmax: float
):
    pp_pt_loc, pp_z_loc = prefix, prefix
    aa_pt_loc, aa_z_loc = prefix, prefix
    if run_type == "with formation time":
        pp_pt_loc += "pp_2p76/jet_FF_pT.csv"
        pp_z_loc += "pp_2p76/jet_FF_z.csv"
        aa_pt_loc += "final_runs/rset_{r}/2p76/{c}/jet_FF_pT.csv"
        aa_z_loc += "final_runs/rset_{r}/2p76/{c}/jet_FF_z.csv"
    else:
        pp_pt_loc += "v2/production_v2/martini_results/pp/jet_FF_pT.csv"
        pp_z_loc += "v2/production_v2/martini_results/pp/jet_FF_z.csv"
        aa_pt_loc += "v2/production_v2/martini_results/final_PbPb_2p76/rset_{r}/cent_{c}/jet_FF_pT.csv"
        aa_z_loc += "v2/production_v2/martini_results/final_PbPb_2p76/rset_{r}/cent_{c}/jet_FF_z.csv"

    pp_pT, pp_z = get_pp_jet_ff(pp_pt_loc, pp_z_loc, pTmin, pTmax, zmin, zmax)
    aa_pT, aa_z = {}, {}
    for r in [1, 2, 3]:
        aa_pT[r] = {}
        aa_z[r] = {}
        for cent in ["0_5", "5_10"]:
            tmp_pt, tmp_z = get_aa_jet_ff(
                aa_pt_loc.format(r=r, c=cent),
                aa_z_loc.format(r=r, c=cent),
                pTmin=pTmin,
                pTmax=pTmax,
                zmin=zmin,
                zmax=zmax,
            )
            aa_pT[r][cent] = tmp_pt
            aa_z[r][cent] = tmp_z

        ## combine the centralities of the aa systems
        tmp_pt = 0.5 * (aa_pT[r]["0_5"]["FF(pT)"] + aa_pT[r]["5_10"]["FF(pT)"])
        tmp_z = 0.5 * (aa_z[r]["0_5"]["FF(z)"] + aa_z[r]["5_10"]["FF(z)"])
        dtmp_pt = 0.5 * np.sqrt(
            aa_pT[r]["0_5"]["dFF(pT)"] ** 2 + aa_pT[r]["5_10"]["dFF(pT)"] ** 2
        )
        dtmp_z = 0.5 * np.sqrt(
            aa_z[r]["0_5"]["dFF(z)"] ** 2 + aa_z[r]["5_10"]["dFF(z)"] ** 2
        )
        aa_pT[r]["0_10"] = pd.DataFrame(
            {
                "pT": aa_pT[r]["0_5"]["pT"],
                "dpT": aa_pT[r]["0_5"]["dpT"],
                "FF(pT)": tmp_pt,
                "dFF(pT)": dtmp_pt,
            }
        )
        aa_z[r]["0_10"] = pd.DataFrame(
            {
                "z": aa_z[r]["0_5"]["z"],
                "dz": aa_z[r]["0_5"]["dz"],
                "FF(z)": tmp_z,
                "dFF(z)": dtmp_z,
            }
        )

    ## Now compute the fragmentation function ratios
    ff_ratio_pt, ff_ratio_z = {}, {}
    for r in aa_z:
        spec_aa_pt, spec_aa_z = aa_pT[r]["0_10"], aa_z[r]["0_10"]

        ratio_pt = spec_aa_pt["FF(pT)"] / pp_pT["FF(pT)"]
        dratio_pt = ratio_pt * np.sqrt(
            (spec_aa_pt["dFF(pT)"] / spec_aa_pt["FF(pT)"]) ** 2
            + (pp_pT["dFF(pT)"] / pp_pT["FF(pT)"]) ** 2
        )

        ratio_z = spec_aa_z["FF(z)"] / pp_z["FF(z)"]
        dratio_z = ratio_z * np.sqrt(
            (spec_aa_z["dFF(z)"] / spec_aa_z["FF(z)"]) ** 2
            + (pp_z["dFF(z)"] / pp_z["FF(z)"]) ** 2
        )
        ff_ratio_pt[r] = pd.DataFrame(
            {
                "x": pp_pT["pT"],
                "dx": pp_pT["dpT"] / 2,
                "ratio": ratio_pt,
                "dratio": dratio_pt,
            }
        )
        ff_ratio_z[r] = pd.DataFrame(
            {
                "x": pp_z["z"],
                "dx": pp_z["dz"] / 2,
                "ratio": ratio_z,
                "dratio": dratio_z,
            }
        )
    return ff_ratio_pt, ff_ratio_z
