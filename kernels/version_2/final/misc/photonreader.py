from numpy import pi, sqrt
from pandas import read_csv, DataFrame
from typing import List

class JetMediumReader:
    def __init__(self, **kw_args):
        """
        Inputs:
            - centralities
            - rate_sets
            - fname_template
            - oversampling
        """
        self.xsec = kw_args["xsec"]
        self.centralities = kw_args["centralities"]
        self.rate_sets = kw_args["rate_sets"]
        self.fname_tmp = kw_args["fname_template"]
        self.oversampling = kw_args["oversampling"]
        self.mult = kw_args["mult"]
        self.spectra = {}
        self._populate_spectra()
        self._construct_0_to_20()
        ## spectra is now full and has the 0-20 % centrality

    def _populate_spectra(self):
        for rate in self.rate_sets:
            self.spectra[rate] = {}
            for c in self.centralities:
                fname = self.fname_tmp.format(r=rate, c=c)
                pTmin, pTmax = (0, 19)
                tmp = read_csv(fname, comment="#")
                if "pTmin" in tmp:
                    tmp["pT"]  = 0.5 * (tmp["pTmin"] + tmp["pTmax"])
                    tmp["dpT"] = tmp["pTmax"] - tmp["pTmin"]
                else:
                    tmp["pT"]  = 0.5 * (tmp["ptmin"] + tmp["ptmax"])
                    tmp["dpT"] = tmp["ptmax"] - tmp["ptmin"]
                tmp = tmp[tmp["pT"].between(pTmin, pTmax)]
                for col in ["conv", "dconv", "brem", "dbrem"]:
                    tmp[col] /= self.xsec * self.oversampling
                    tmp[col] *= self.mult[c]
                tmp["conv"]  /= 2 * pi * 2 * 0.8 * tmp["pT"] * tmp["dpT"]
                tmp["dconv"] /= 2 * pi * 2 * 0.8 * tmp["pT"] * tmp["dpT"]
                tmp["brem"]  /= 2 * pi * 2 * 0.8 * tmp["pT"] * tmp["dpT"]
                tmp["dbrem"] /= 2 * pi * 2 * 0.8 * tmp["pT"] * tmp["dpT"]
                tmp["total"]  = tmp["conv"] + tmp["brem"]
                tmp["dtotal"] = sqrt(
                    tmp["dconv"] * tmp["dconv"] + tmp["dbrem"] * tmp["dbrem"]
                )
                self.spectra[rate][c] = tmp

    def _construct_0_to_20(self):
        for rate in self.rate_sets:
            spec = self.spectra[rate]
            tmp = {}
            tmp["pT"] = spec["0_5"]["pT"]
            tmp["conv"] = 0.25 * (
                spec["0_5"]["conv"] + spec["5_10"]["conv"] + spec["10_20"]["conv"]
            )
            tmp["dconv"] = 0.25 * sqrt(
                spec["0_5"]["dconv"] ** 2
                + spec["5_10"]["dconv"] ** 2
                + spec["10_20"]["dconv"] ** 2
            )

            tmp["brem"] = 0.25 * (
                spec["0_5"]["brem"] + spec["5_10"]["brem"] + spec["10_20"]["brem"]
            )
            tmp["dbrem"] = 0.25 * sqrt(
                spec["0_5"]["dbrem"] ** 2
                + spec["5_10"]["dbrem"] ** 2
                + spec["10_20"]["dbrem"] ** 2
            )

            tmp["total"] = 0.25 * (
                spec["0_5"]["total"] + spec["5_10"]["total"] + spec["10_20"]["total"]
            )
            tmp["dtotal"] = 0.25 * sqrt(
                spec["0_5"]["dtotal"] ** 2
                + spec["5_10"]["dtotal"] ** 2
                + spec["10_20"]["dtotal"] ** 2
            )
            self.spectra[rate]["0_20"] = DataFrame(tmp)

    def get_spectra(self):
        return self.spectra


#class NonJetMediumReader:
#    """
#        Class to read in the non jet-medium photons
#            * reads in the photons channels that are non-jetmedium
#            * harmonizes the x axes across all of them to the provided
#            axis
#    """
#    def __init__(self, xvalues: List[float]):
#        self.xvals = xvalues
#        self.spectra = {}
#    
#    def read_photon_channel(fname : str, channel_name : str):
#        


class JETSCAPEReader:
    def __init__(self, **kw_args):
        """
        Inputs:
            - centralities
            - modelname 
            - fname_template
            - oversampling factors
            - fname_template
            - mult
        """
        self.xsec = kw_args["xsec"]
        self.centralities = kw_args["centralities"]
        self.modelname = kw_args["modelname"]
        self.fname_tmp = kw_args["fname_template"]
        self.oversampling = kw_args["oversampling"]
        self.mult = kw_args["mult"]
        self.spectra = {}
        self._populate_spectra()
        self._construct_0_to_20()
        ## spectra is now full and has the 0-20 % centrality

    def _populate_spectra(self):
        for cent, oversample in zip(self.centralities, self.oversampling):
            fname = self.fname_tmp.format(model=self.modelname, c=cent)
            pTmin, pTmax = (0, 19)
            tmp = read_csv(fname, comment="#")
            if "pTmin" in tmp:
                tmp["pT"]  = 0.5 * (tmp["pTmin"] + tmp["pTmax"])
                tmp["dpT"] = tmp["pTmax"] - tmp["pTmin"]
            else:
                tmp["pT"]  = 0.5 * (tmp["ptmin"] + tmp["ptmax"])
                tmp["dpT"] = tmp["ptmax"] - tmp["ptmin"]
            tmp = tmp[tmp["pT"].between(pTmin, pTmax)]
            for col in ["conv", "dconv", "brem", "dbrem"]:
                tmp[col] /= self.xsec * oversample
                tmp[col] *= self.mult[cent]
            tmp["conv"]  /= 2 * pi * 2 * 0.8 * tmp["pT"] * tmp["dpT"]
            tmp["dconv"] /= 2 * pi * 2 * 0.8 * tmp["pT"] * tmp["dpT"]
            tmp["brem"]  /= 2 * pi * 2 * 0.8 * tmp["pT"] * tmp["dpT"]
            tmp["dbrem"] /= 2 * pi * 2 * 0.8 * tmp["pT"] * tmp["dpT"]
            tmp["total"]  = tmp["conv"] + tmp["brem"]
            tmp["dtotal"] = sqrt(
                tmp["dconv"] * tmp["dconv"] + tmp["dbrem"] * tmp["dbrem"]
            )
            self.spectra[cent] = tmp

    def _construct_0_to_20(self):
        spec = self.spectra
        tmp = {}
        tmp["pT"] = spec["00-05"]["pT"]
        tmp["conv"] = 0.25 * (
            spec["00-05"]["conv"] + spec["05-10"]["conv"] + spec["10-20"]["conv"]
        )
        tmp["dconv"] = 0.25 * sqrt(
              spec["00-05"]["dconv"] ** 2
            + spec["05-10"]["dconv"] ** 2
            + spec["10-20"]["dconv"] ** 2
        )

        tmp["brem"] = 0.25 * (
            spec["00-05"]["brem"] + spec["05-10"]["brem"] + spec["10-20"]["brem"]
        )
        tmp["dbrem"] = 0.25 * sqrt(
              spec["00-05"]["dbrem"] ** 2
            + spec["05-10"]["dbrem"] ** 2
            + spec["10-20"]["dbrem"] ** 2
        )

        tmp["total"] = 0.25 * (
            spec["00-05"]["total"] + spec["05-10"]["total"] + spec["10-20"]["total"]
        )
        tmp["dtotal"] = 0.25 * sqrt(
              spec["00-05"]["dtotal"] ** 2
            + spec["05-10"]["dtotal"] ** 2
            + spec["10-20"]["dtotal"] ** 2
        )
        self.spectra["0_20"] = DataFrame(tmp)

    def get_spectra(self):
        return self.spectra