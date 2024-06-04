#include <iostream>
#include <fstream>
#include <iomanip>
#include <chrono>
#include <algorithm>
#include <cmath>
#include "analysis.h"
int main(int argc, char* argv[]){

    string pT_str = argv[1];
    double pT = stod(pT_str);
    string nevents     = argv[2];
    string save_loc    = argv[3];
    string setup_fname = argv[4];

    int numEvents = std::stoi(nevents);
    
    // Prepare MARTINI
    MARTINI martini;
    martini.readFile(setup_fname);
    martini.readString("General:JetPTMin = " + pT_str);
    martini.readString("General:JetPTMax = " + pT_str);
    martini.init(0);
    martini.settings.listChanged(); // list the changed parameters. 
    vector<Parton> * plist = new vector<Parton>; //pointer to the parton list vector
    vector<Source> * slist = NULL; // pointer to the source list, null for now

    bool medium_evolution = martini.returnEvolution();
    int mt; // maximal time steps
    double maxTime = martini.returnMaxTime();
    double dtfm = martini.returnDtfm();// dt in femtometers 
    mt = static_cast<int>(maxTime/dtfm+0.0001); //max number of steps
    // cout << maxTime << ", "<<dtfm << ", " << mt << endl;

    float first = 0.05, last = 1.05, step = 0.05;
    // z is the ratio of the pT of the parton to the 
    // input pT. So we're measuring with respect to
    // the "leading" parton.
    int num_z_bins = ceil((last-first)/step)+2;
    vector<int> chosen_times = {4, 8, 12, 16, 20, 24, 28, 32, 
                               36, 40, 44, 48, 52, 56, 60, 64, 68, 
                               72, 76, 80, 84, 88, 92, 96, 100, 104, 
                               108, 112, 116, 120, 124, 128, 132, 136, 140};
    
    int num_time_hists = chosen_times.size();


    vector<double> zbins(num_z_bins, 0.0);
    for (size_t iz = 0; iz < num_z_bins; iz++)
        zbins[iz] = iz * step;

    vector<vector<double>> quarks;
    vector<vector<double>> gluons;
    quarks.reserve(num_time_hists);
    gluons.reserve(num_time_hists);

    for (size_t i = 0; i < num_time_hists; i++){
        vector<double> tmp(num_z_bins, 0.0);
        quarks.push_back(tmp);

        vector<double> tmp2(num_z_bins, 0.0);
        gluons.push_back(tmp2);
    }
    
    int counter       = 0;
    int event_counter = 0; 

    auto start_program = std::chrono::steady_clock::now();
    auto end_program   = std::chrono::steady_clock::now();
    auto elapsed_program = std::chrono::duration_cast<std::chrono::minutes>(end_program - start_program);
    int time_limit = 5*55; // minutes
    int time_count_prog;
    int prop_time = 0;
    double tau = 0;
    while( event_counter < numEvents ) {
        // put the time check both here and at the end
        elapsed_program = std::chrono::duration_cast<std::chrono::minutes>(end_program - start_program);
        time_count_prog = elapsed_program.count();
        if ( time_count_prog > time_limit - 5 ) { // give it 5 minutes to wrap up
            numEvents = event_counter; // we only had this many events        
            cout << "Have run for "<< time_count_prog << " minutes. ";
            cout << "Exiting the event loop."<<endl;
            break;
        }
        plist->clear();
        martini.generateEvent(plist);
        counter = 0;
        // cout << "num_time_hists: " << num_time_hists << endl;
        if (medium_evolution){// evolve in medium if settings want that
            for(int i=0; i<mt; i++){ // loop over all time steps 
                counter = martini.evolve(plist, slist, counter, i);
                counter += 1;
                if (i % 200 == 0){
                    tau = 0.4 + i * dtfm;
                    prop_time = (int) (10 * tau);
                    // cout << "prop_time: "<<prop_time<<endl;
                    auto pos = std::find(chosen_times.begin(), chosen_times.end(), prop_time);
                    if ( pos != chosen_times.end() ){
                        // cout << "prop time "<< prop_time << " is in the chosen times! "<<endl;
                        bin_partons_by_z(plist, gluons[pos - chosen_times.begin()], 
                                                quarks[pos - chosen_times.begin()], zbins, pT, tau);
                    }
                }
            }
        }
        event_counter++;
        end_program = std::chrono::steady_clock::now();
    }
    // done with the loop
    martini.pythia.stat();
    double sigmaGen = martini.pythia.info.sigmaGen();

    // Save to file:
    // cout << "Num histograms: " << gluons.size() << endl;
    // cout << "Num bins in a histogram: " << gluons[0].size() << endl;
    // cout << "num_time_hists: "<<num_time_hists << endl;
    // cout << "num_z_bins: "<< num_z_bins << endl;
    double n, dn;
    string fname;
    for (int i = 0; i < num_time_hists; i++){
        fname = save_loc + "tau_" + to_string(i) + "_partons_gun_" + pT_str + ".csv";
        fstream out( fname.c_str(), ios::out);
        out << "#nevents " << numEvents << " sigmaGen ";
        out << sigmaGen << " eta "<< PARTON_ETA_CUT <<endl;
        out << "#tau " << setprecision(3) << chosen_times[i]/10.0 << endl;
        out <<"zmin,zmax,Nglue,Nquark"<<endl;
        for (size_t ibin = 0; ibin < num_z_bins-1; ibin++){
            out << zbins[ibin] << ",";        
            out << zbins[ibin+1];
            n  = gluons[i][ibin];
            n = n > 0 ? n : 0;
            out << "," << n;
            n = quarks[i][ibin];
            n = n > 0 ? n : 0;
            out << ","<<n;
            out << endl;
        }
    }
    cout << "Have run for "<< time_count_prog << " minutes."<<endl;
    cout << "sigmaGen : "  << sigmaGen << endl;
    cout << "numEvents : " << numEvents <<endl;
    return 0;
}
