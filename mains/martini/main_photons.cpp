#include <iostream>
#include <fstream>
#include <iomanip>
#include <chrono>
#include "analysis.h"
int main(int argc, char* argv[]){
    string subrun      = argv[1];
    string pTmin_str   = argv[2];
    string pTmax_str   = argv[3];
    string nevents     = argv[4];
    string fit_run     = argv[5];
    string save_loc    = argv[6];
    string setup_fname = argv[7];
    // cout << subrun<< " ** "<<nevents << " ** "<<save_loc<<" ** " << setup_fname<<endl;
    // number of events, jet pT info
    int numEvents = std::stoi(nevents);
    int fit_run_int = std::stoi(fit_run);
    bool is_fit_run = fit_run_int == 1 ? true : false; 
    
    // Prepare MARTINI
    MARTINI martini;
    martini.readFile(setup_fname);
    cout << "subrun #: " << subrun<<endl;
    martini.readString("General:JetPTMin = " + pTmin_str);
    martini.readString("General:JetPTMax = " + pTmax_str);
    martini.pythia.particleData.readString("111:mayDecay = off");
    martini.init(0);
    martini.settings.listChanged(); // list the changed parameters. 
    vector<Parton> * plist = new vector<Parton>; //pointer to the parton list vector
    vector<Source> * slist = NULL; // pointer to the source list, null for now

    bool medium_evolution = false; //martini.returnEvolution();
    int produce_photons = 1; //martini.returnPhotonSwitch();    
    int do_fragmentation = 0;//martini.returnFragmentationSwitch();
    int mt; // maximal time steps
    double maxTime = martini.returnMaxTime();
    double dtfm = martini.returnDtfm();// dt in femtometers 
    mt = static_cast<int>(maxTime/dtfm+0.0001); //max number of steps

    // PHOTON SPECTRA HISTOGRAM
    vector<double> photon_spec_bins  = {0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5,
                                        4.0, 4.5, 5.0, 5.5, 6.0, 6.7, 7.0, 7.5,
                                        8.0, 8.5, 9.0, 9.5, 10., 12., 14., 16.,
                                        18., 20., 22., 24., 26., 28., 30., 35., 40., 45., 50.,
                                        55., 60., 65., 70., 75., 80., 85., 90., 100., 999};
    size_t nbins_photons = photon_spec_bins.size()-1; 
    vector<double> prompt_hist(nbins_photons, 0.0);
    vector<double> conv_hist(nbins_photons, 0.0);
    vector<double> brem_hist(nbins_photons, 0.0);
    int counter       = 0;
    int event_counter = 0; 
    int hadronized    = 0;

    auto start_program = std::chrono::steady_clock::now();
    auto end_program   = std::chrono::steady_clock::now();
    auto elapsed_program = std::chrono::duration_cast<std::chrono::minutes>(end_program - start_program);
    int time_limit = 3*55; // minutes
    int time_count_prog;
    Event hardEvent, recoilEvent, holeEvent;

    while( event_counter < numEvents ){
        // put the time check both here and at the end
        elapsed_program = std::chrono::duration_cast<std::chrono::minutes>(end_program - start_program);
        time_count_prog = elapsed_program.count();
        if ( time_count_prog > time_limit - 5 ) { // give it 5 minutes to wrap up
            numEvents = event_counter; // we only had this many events        
            cout << "Have run for "<< time_count_prog << " minutes. ";
            cout <<"Exiting the event loop."<<endl;
            break;
        }
        plist->clear();
        martini.generateEvent(plist);
        counter = 0;
        if (medium_evolution){// evolve in medium if settings want that
            for(int i=0; i<mt; i++){ // loop over all time steps 
                counter = martini.evolve(plist, slist, counter, i);
                counter+=1;
            }
        }
        if (produce_photons and not is_fit_run)
            bin_photons(plist, prompt_hist, conv_hist, brem_hist, photon_spec_bins);

        event_counter++;
        hadronized += 1;
        end_program = std::chrono::steady_clock::now();
    }
    // done with the loop
    martini.pythia.stat();
    double sigmaGen = martini.pythia.info.sigmaGen();

    // Save to file:
    double n, dn;
    // PHOTONS
    if (produce_photons){
        string fname = save_loc + "photons_" + subrun + ".csv";
        fstream output( fname.c_str() ,ios::out);
        output << "#hadronized " << hadronized << " sigmaGen "<< sigmaGen << " eta "<< PHOTON_ETA_CUT<<endl;
        output <<"pTmin,pTmax,prompt,dprompt,conv,dconv,brem,dbrem"<<endl;
        for (size_t ibin = 0; ibin < nbins_photons; ibin++){
            output<< photon_spec_bins[ibin] << ",";
            output<< photon_spec_bins[ibin+1]<<",";
            n  = prompt_hist[ibin];
            dn = sqrt(n);
            output << n << ",";
            output << dn << ",";
            n  = conv_hist[ibin];
            dn = sqrt(n);
            output << n << ",";
            output << dn << ",";
            n  = brem_hist[ibin];
            dn = sqrt(n);
            output << n << ",";
            output << dn << endl;
        }
        output.close();
    }
    cout << "Have run for "<< time_count_prog << " minutes."<<endl;
    cout << "sigmaGen : "  << sigmaGen << endl;
    cout << "numEvents : " << numEvents <<endl;
    cout << "Num hadronized: "  << hadronized  << " --> " << hadronized/(double)numEvents << endl;
    return 0;
}
