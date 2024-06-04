#include <iostream>
#include <fstream>
#include <iomanip>
#include <chrono>
#include "analysis.h"
#define JET_MIN_PT 20
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
    martini.init(0);
    martini.settings.listChanged(); // list the changed parameters. 
    vector<Parton> * plist = new vector<Parton>; //pointer to the parton list vector
    vector<Source> * slist = NULL; // pointer to the source list, null for now

    bool medium_evolution = martini.returnEvolution();
    int produce_photons = martini.returnPhotonSwitch();    
    int do_fragmentation = martini.returnFragmentationSwitch();
    int mt; // maximal time steps
    double maxTime = martini.returnMaxTime();
    double dtfm = martini.returnDtfm();// dt in femtometers 
    mt = static_cast<int>(maxTime/dtfm+0.0001); //max number of steps

    // CHARGED HADRON HISTOGRAM
    vector<double> charged_spec_bins = {1 , 2  , 3  , 4  , 5  , 6  , 7  , 8  , 9  ,
                                 10 , 12 , 14 , 16 , 18 , 20 , 23 , 26 , 30 ,
                                 35 , 40 , 45 , 50 , 55 , 60 , 70 , 80 , 90 ,
                                 100, 120, 140, 160, 180, 200, 250, 300, 350, 
                                 400, 500, 9999};
    size_t nbins_ch_hads = charged_spec_bins.size()-1; 
    vector<double> charged_hist(nbins_ch_hads, 0.0);
    // Parton Spectra will have the same binning as charged hadrons.
    // Histograms for fermions and gluons one for events that failed at 
    // hadronization and one for those that did not.

    // JET SPECTRA HISTOGRAM
    vector<double> jet_spec_bins = {20 , 22, 24, 26, 29, 34, 40, 47, 55, 64 ,
                                  74 ,84 ,97 ,114,133,153,174,196,220,245 ,272,
                                  300,350,400,500,550,600,700,800,900,1000, 9999};
    size_t nbins_jets    = jet_spec_bins.size()-1;
    vector<double> jet_radii = {0.2, 0.3, 0.4};
    size_t num_jet_radii = jet_radii.size();

    vector<vector<double>> jet_hists;
    vector<JetDefinition> jet_defs;

    jet_hists.reserve(num_jet_radii);
    jet_defs.reserve(num_jet_radii);

    vector<vector<double>> partonic_jets_success;
    partonic_jets_success.reserve(num_jet_radii);
    
    for (auto r: jet_radii){
        vector<double> temp(nbins_jets, 0.0);
        jet_hists.push_back(temp);
        
        vector<double> tmp_part_succ(nbins_jets, 0.0);
        partonic_jets_success.push_back(tmp_part_succ);

        JetDefinition jet_def(antikt_algorithm, r);
        jet_defs.push_back(jet_def);
    }
    // PHOTON SPECTRA HISTOGRAM
    vector<double> photon_spec_bins  = { 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 
                                        14, 15, 17, 19, 21, 23, 25, 30, 40, 45, 9999};
    size_t nbins_photons = photon_spec_bins.size()-1; 
    vector<double> prompt_hist(nbins_photons, 0.0);
    vector<double> conv_hist(nbins_photons, 0.0);
    vector<double> brem_hist(nbins_photons, 0.0);

    // FRAGMENTATION FUNCTION HISTOGRAMS
    vector<double> zbins = {0.0, 0.01, 0.016, 0.025, 0.04, 0.063, 0.1, 0.16, 
                            0.25, 0.4, 0.63, 1.0};
    vector<double> FFpTbins = {0.0, 1.0, 1.6, 2.5, 4, 6.3, 10, 16, 25, 40, 63, 100, 9999};
    size_t nbins_z = zbins.size()-1;
    size_t nbins_FFpTbins = FFpTbins.size()-1;
    vector<double> FFz(nbins_z, 0.0);
    vector<double> FFpT(nbins_FFpTbins, 0.0);
    // SHAPE HISTOGRAM INFO
    vector<double> rbins = {0, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45};
    size_t num_rbins = rbins.size()-1;
    vector<double> shape(num_rbins, 0.0);
    // Also keep jet shape information for partonic jets,
    // only for R=0.4 at mid rapidity
    vector<double> shape_gluonic_success(num_rbins, 0.0);
    vector<double> shape_fermionic_success(num_rbins, 0.0);
    
    int counter       = 0;
    int event_counter = 0; 
    int hadronized    = 0;
    int num_processed_jets_FF = 0;
    int num_processed_jets_shape = 0;
    int num_processed_success_partonic_shape = 0;
    int tmp_FF=0, tmp_shape=0, tmp_shape_part=0;

    auto start_program = std::chrono::steady_clock::now();
    auto end_program   = std::chrono::steady_clock::now();
    auto elapsed_program = std::chrono::duration_cast<std::chrono::minutes>(end_program - start_program);
    int time_limit = 5*55; // minutes
    int time_count_prog;
    Event hardEvent, holeEvent;
    //recoilEvent ;

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

        hardEvent.clear();
        //recoilEvent.clear();
        holeEvent.clear();
        martini.pythia.event.clear();

        // Attempt to fragment the hard particles
        if (not martini.fragmentation(plist, HARD))
            continue;
        hardEvent = martini.pythia.event;
        if (not martini.fragmentation(plist, HOLE))
            continue;
        holeEvent = martini.pythia.event;

        if (not is_fit_run){
            vector<PseudoJet> success_partonic_pseudojets;
            vector<PseudoJet> success_partonic_holes;
            cluster_prep_partonic(plist, success_partonic_pseudojets, HARD);
            cluster_prep_partonic(plist, success_partonic_holes, HOLE);
            // cluster 
            for (size_t ir = 0; ir < num_jet_radii; ir++){
                ClusterSequence cs(success_partonic_pseudojets, jet_defs[ir]);
                vector<PseudoJet> clustered_jets = sorted_by_pt(cs.inclusive_jets(JET_MIN_PT));
                subtract_holes_from_jets(clustered_jets, success_partonic_holes, jet_radii[ir]);
                bin_jets(clustered_jets, partonic_jets_success[ir], jet_spec_bins);
                if (ir == 2){
                    tmp_shape_part = bin_partonic_shape(clustered_jets, shape_gluonic_success, rbins, GLUON);
                    bin_partonic_shape(clustered_jets, shape_fermionic_success, rbins, QUARK);
                    num_processed_success_partonic_shape += tmp_shape_part;
                }
            }
        }
        // everything has been successful. Compute observables and bin it all.
        // bin charged hadrons
        bin_chgd_hads(hardEvent, charged_hist, charged_spec_bins, 1);
        bin_chgd_hads(holeEvent, charged_hist, charged_spec_bins, -1);

        vector<PseudoJet> forFJ;
        vector<PseudoJet> holes;

        cluster_prep_inclusive(hardEvent, forFJ);
        cluster_prep_inclusive(holeEvent, holes);//subtracted after jets are clustered

        //Now loop over the jet radii
        for (size_t ir = 0; ir < num_jet_radii; ir++){
            // For this jet definition, cluster the jets 
            ClusterSequence cs(forFJ, jet_defs[ir]);
            vector<PseudoJet> clustered_jets = sorted_by_pt(cs.inclusive_jets(JET_MIN_PT));
            subtract_holes_from_jets(clustered_jets, holes, jet_radii[ir]);
            bin_jets(clustered_jets, jet_hists[ir], jet_spec_bins);
            if (ir == 2 and not is_fit_run){
                // Jet Fragmentation Functions are computed only for R = 0.4
                tmp_FF = bin_jet_fragmentation(clustered_jets, FFz, zbins, FFpT, FFpTbins);
                num_processed_jets_FF += tmp_FF;
            }
            if (ir == 1){
                vector<PseudoJet> forFJ_wcut;
                cluster_prep_inclusive_wth_trk_cut(hardEvent, forFJ_wcut, trk_cut_shape);
                ClusterSequence tmp(forFJ_wcut, jet_defs[ir]);
                vector<PseudoJet> clustered_shapes = sorted_by_pt(tmp.inclusive_jets(JET_MIN_PT));
                subtract_holes_from_jets(clustered_shapes, holes, jet_radii[ir]);
                tmp_shape = bin_jet_shape(clustered_shapes, shape, rbins);
                num_processed_jets_shape += tmp_shape;
            }
        }
        hadronized += 1;
        end_program = std::chrono::steady_clock::now();
    }
    // done with the loop
    martini.pythia.stat();
    double sigmaGen = martini.pythia.info.sigmaGen();

    // Save to file:
    double n, dn;
    if (do_fragmentation){
        // CHARGED HADRONS
        string fname_charged_hadrons = save_loc + "charged_hadrons_" + subrun + ".csv";
        fstream charged_hadron_file( fname_charged_hadrons.c_str() ,ios::out);
        charged_hadron_file << "#hadronized " << hadronized;
        charged_hadron_file << " sigmaGen "<< sigmaGen << " eta "<< CH_HAD_ETA_CUT<<endl;
        charged_hadron_file << "pTmin,pTmax,N,dN"<<endl;
        for (size_t ibin = 0; ibin < nbins_ch_hads; ibin++){
            charged_hadron_file << charged_spec_bins[ibin] << ",";
            charged_hadron_file << charged_spec_bins[ibin+1] << ",";
            n  = charged_hist[ibin];
            n = n>0 ? n : 0;
            dn = sqrt(n);
            charged_hadron_file << n << ",";               
            charged_hadron_file << dn << endl;
        }
        charged_hadron_file.close();
        // HADRONIC JETS
        string fname_jets = save_loc + "jets_" + subrun + ".csv";
        fstream jets_file( fname_jets.c_str() ,ios::out);
        jets_file << "#hadronized " << hadronized << " sigmaGen ";
        jets_file << sigmaGen << " eta "<< JET_ETA_CUT <<endl;
        jets_file <<"pTmin,pTmax,N_0p2,dN_0p2,N_0p3,dN_0p3,N_0p4,dN_0p4"<<endl;
        for (size_t ibin = 0; ibin < nbins_jets; ibin++){
            jets_file << jet_spec_bins[ibin] << ",";        
            jets_file << jet_spec_bins[ibin+1];
            for (size_t ir = 0; ir < num_jet_radii; ir++){
                n  = jet_hists[ir][ibin];
                n = n > 0 ? n : 0;
                dn = sqrt(n);
                jets_file<< "," << n;
                jets_file<< "," << dn;
            }
            jets_file<<endl;
        }
        jets_file.close();
        if (not is_fit_run){
            // JET Fragmentation Function: FF(z)
            string fname_FFz = save_loc + "jet_FF_z_" + subrun + ".csv";
            fstream FFz_out(fname_FFz.c_str(), ios::out);
            FFz_out << "#hadronized " << hadronized<<" sigmaGen ";
            FFz_out << sigmaGen << " num_jets "<< num_processed_jets_FF;
            FFz_out << " JET_PTMin " << JET_PTMIN_FF << " JET_PTMax ";
            FFz_out << JET_PTMAX_FF;
            FFz_out << " Jet_R " << "0.4" << endl;
            FFz_out << "zmin,zmax,N,dN"<<endl;
            for (size_t ibin = 0; ibin < nbins_z; ibin++){
                FFz_out << zbins[ibin] << ",";
                FFz_out << zbins[ibin+1] << ",";
                n = FFz[ibin];
                dn = sqrt(n);
                FFz_out << n << ",";
                FFz_out << dn<< endl;
            }
            FFz_out.close();
            // JET Fragmentation Function: FF(pT)
            string fname_FFpT = save_loc + "jet_FF_pT_" + subrun + ".csv";
            fstream FFpT_out(fname_FFpT.c_str(), ios::out);
            FFpT_out << "#hadronized " << hadronized<<" sigmaGen ";
            FFpT_out << sigmaGen << " num_jets "<< num_processed_jets_FF;
            FFpT_out << " JET_PTMin " << JET_PTMIN_FF << " JET_PTMax ";
            FFpT_out << JET_PTMAX_FF;
            FFpT_out << " Jet_R " << "0.4" << endl;
            FFpT_out << "pTmin,pTmax,N,dN" << endl;
            for (size_t ibin = 0; ibin < nbins_FFpTbins; ibin++){
                FFpT_out << FFpTbins[ibin] << ",";
                FFpT_out << FFpTbins[ibin+1] << ",";
                n = FFpT[ibin];
                dn = sqrt(n);
                FFpT_out << n << ",";
                FFpT_out << dn<< endl;
            }
            FFpT_out.close();
        }
        // JET Shape 
        string fname_shape = save_loc + "jet_shape_" + subrun + ".csv";
        fstream fshape(fname_shape.c_str(), ios::out);
        fshape << "#hadronized " << hadronized<<" sigmaGen ";
        fshape << sigmaGen << " num_jets "<< num_processed_jets_shape;
        fshape << " JET_PTMin " << JET_PTMIN_SHAPE;
        fshape << " Jet_R " << "0.3 withcut 1.0" << endl;
        fshape << "rmin,rmax,N,dN"<<endl;
        for (size_t ibin = 0; ibin < num_rbins; ibin++){
            fshape << rbins[ibin] << ",";
            fshape << rbins[ibin+1] << ",";
            n = shape[ibin];
            dn = sqrt(n);
            fshape << n << ",";
            fshape << dn<< endl;
        }
        fshape.close();
        // write the parton spectra (quarks and gluons) of failed and 
        // successful 
        if (not is_fit_run){
            // write partonic jets
            string fname_jets_success = save_loc + "successful_parton_jets_" + subrun + ".csv";
            fstream jets_file_success( fname_jets_success.c_str(), ios::out);
            jets_file_success << "#hadronized " << hadronized << " sigmaGen ";
            jets_file_success << sigmaGen << " eta "<< JET_ETA_CUT <<endl;
            jets_file_success <<"pTmin,pTmax,N_0p2,dN_0p2,N_0p3,dN_0p3,N_0p4,dN_0p4"<<endl;
            for (size_t ibin = 0; ibin < nbins_jets; ibin++){
                jets_file_success << jet_spec_bins[ibin] << ",";        
                jets_file_success << jet_spec_bins[ibin+1];
                for (size_t ir = 0; ir < num_jet_radii; ir++){
                    n = partonic_jets_success[ir][ibin];
                    n = n > 0 ? n : 0;
                    dn = sqrt(n);
                    jets_file_success<< "," << n;
                    jets_file_success<< "," << dn;
                }
                jets_file_success<<endl;
            }
            jets_file_success.close();

            // Partonic Jet Shape:
            // successful partonic shapes
            string fname_shape_part_success = save_loc + "success_partonic_jet_shape__" + subrun + ".csv";
            fstream fshape_ps(fname_shape_part_success.c_str(), ios::out);
            fshape_ps << "#hadronized " << hadronized <<" sigmaGen ";
            fshape_ps << sigmaGen << " num_jets "<< num_processed_success_partonic_shape;
            fshape_ps << " JET_PTMin " << JET_PTMIN_SHAPE;
            fshape_ps << " Jet_R " << "0.4 no cut" << endl;
            fshape_ps << "rmin,rmax,fermion,dfermion,gluon,dgluon"<<endl;
            for (size_t ibin = 0; ibin < num_rbins; ibin++){
                fshape_ps << rbins[ibin] << ",";
                fshape_ps << rbins[ibin+1] << ",";
                n = shape_fermionic_success[ibin];
                dn = sqrt(n);
                fshape_ps << n << ",";
                fshape_ps << dn<< ",";
                n = shape_gluonic_success[ibin];
                dn = sqrt(n);
                fshape_ps << n << ",";
                fshape_ps << dn<< endl;
            }
            fshape_ps.close();
        }
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
    }
    cout << "Have run for "<< time_count_prog << " minutes."<<endl;
    cout << "sigmaGen : "  << sigmaGen << endl;
    cout << "numEvents : " << numEvents <<endl;
    cout << "Num hadronized: "  << hadronized  << " --> " << hadronized/(double)numEvents << endl;
    return 0;
}
