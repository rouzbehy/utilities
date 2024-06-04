from numpy import array, sqrt, tanh, cos, pi
import os

def get_xsection(loc):
    res = -1
    if os.path.exists(loc + "sigmaGen.dat"):
        ## sigmaGen.dat is there, simply read it
        with open(loc + "sigmaGen.dat", 'r') as f:
            line = f.readline()
            (sigmaGen, _) = line.split('\t')
            res = float(sigmaGen)
    else:
        fname = loc + "evolution_result.dat"
        if not os.path.exists(fname):
            print(f"this {loc} run does not even an evolution file. No Sigma Gen")
        else:
            for line in open(fname, 'r'):
                if 'sigmaGen' in line:
                    line = line.split(' ')
                    sigmaGen = line[2]
            ## finally the last sigmaGen is the most accurate one
            res = float(sigmaGen)
    return res
class Histogram:
    def __init__(self, xvals):
       self.bins = array(xvals)
       self.nbins = len(xvals)-1
       self.hist = array([0. for x in range(self.nbins)])
       self.err  = array([0. for x in range(self.nbins)])
       self.nentries = 0

    def fill(self, xval, weight=1):
        for i in range(self.nbins):
            if self.bins[i] <= xval < self.bins[i+1]:
                self.hist[i] += weight
                self.err[i] += weight*weight
                self.nentries += 1
                return True
        return False

    def get_nbins(self):
        return self.nbins

    def multiply_by(self, number):
        self.hist = self.hist*number
        self.err = self.err*(number*number)

    def get_xbins(self):
        return self.bins

    def get_content(self):
        return self.hist, sqrt(self.err)

    def set_content(self, yvals, dyvals):
        self.hist = yvals
        self.err = dyvals*dyvals

    def __getitem__(self, index):
        return self.hist[index], sqrt(self.err[index])

    def __add__(self, other):
        try:
            ohist, oerr = other.get_content()
            new_content = self.hist + ohist
            new_err = sqrt(self.err + oerr*oerr)
        except:
            print(ohist)
            print(self.hist)
            exit(-1)
        newHist = Histogram(self.bins)
        newHist.set_content(new_content, new_err)
        return newHist
            

class Particle:
    def __init__(self, line):
        ## jet id associated with this particle
        self.jetNum = int(line[0].split('\t')[1])
        self.hadNum = int(line[1])## id of this particle within the jet
        self.hadID = int(line[2]) ## PDG particle ID
        self.en  = float(line[3])
        self.pT  = float(line[4]) 
        self.phi = float(line[5]) 
        tmp_eta = line[6]
        if tmp_eta == 'nan':
            self.eta = 100.
        elif tmp_eta == '-nan':
            self.eta = -100.
        else:
            self.eta = float(tmp_eta)
        
class Jet:
    def __init__(self, line):
        self.jetNum = int(line[0].split('\t')[1])
        self.en  = float(line[1]) 
        self.pT  = float(line[2])
        self.phi = float(line[3]) 
        self.eta = float(line[4]) 

    def calc_r(self, hadron):
        delta_phi = abs(self.phi - hadron.phi)
        delta_eta = abs(self.eta - hadron.eta)
        if delta_phi > pi:
            delta_phi = abs(2*pi - delta_phi)
        if delta_phi > pi:
            delta_phi = abs(2*pi - delta_phi)
        return sqrt(delta_phi*delta_phi + delta_eta*delta_eta)

    def calc_z(self, hadron):
        delta_phi = abs(self.phi - hadron.phi)
        delta_eta = abs(self.eta - hadron.eta)
        if delta_phi > pi:
            delta_phi = abs(2*pi - delta_phi)
        if delta_phi > pi:
            delta_phi = abs(2*pi - delta_phi)
        pJet_dot_pJet = (self.pT**2) + (self.en*tanh(self.eta))**2
        pJet_dot_pTrk = self.pT*hadron.pT*cos(delta_phi)+ self.en*hadron.en*tanh(self.eta)*tanh(hadron.eta)
        z = pJet_dot_pTrk / pJet_dot_pJet
        return z
    def calc_del_phi(self, other):
       delta_phi = abs(self.phi - other.phi)
       if delta_phi > pi:
           delta_phi = abs(2*pi - delta_phi)
       if delta_phi > pi:
           delta_phi = abs(2*pi - delta_phi)
       return  delta_phi


if __name__=='__main__':
    xbins = [1, 2, 4, 7]
    hist1 = Histogram(xbins)
    hist2 = Histogram(xbins)
    import random
    for i in range(20):
        hist1.fill(random.random()*12)
        hist2.fill(random.random()*20)
    print(hist1.hist)
    print(hist2.hist)
    hist3 = hist1+ hist2
    print(hist3.hist)

