/*******************************************************************************
 * Copyright (c) The JETSCAPE Collaboration, 2018
 *
 * Modular, task-based framework for simulating all aspects of heavy-ion collisions
 * 
 * For the list of contributors see AUTHORS.
 *
 * Report issues at https://github.com/JETSCAPE/JETSCAPE/issues
 *
 * or via email to bugs.jetscape@gmail.com
 *
 * Distributed under the GNU General Public License 3.0 (GPLv3 or later).
 * See COPYING for details.
 ******************************************************************************/

#include <iostream>
#include <fstream>
#include <sstream>
#include <memory>
#include <chrono>
#include <thread>

#include "gzstream.h"
#include "PartonShower.h"
#include "JetScapeLogger.h"
#include "JetScapeReader.h"
#include "JetScapeBanner.h"
#include "fjcore.hh"

#include <GTL/dfs.h>

using namespace std;
using namespace fjcore;
using namespace Jetscape;

// You could overload here and then simply use ofstream << p;
// ostream & operator<<(ostream & ostr, const fjcore::PseudoJet & jet);
// -------------------------------------

int main(int argc, char** argv)
{
    JetScapeLogger::Instance()->SetDebug(false);
    JetScapeLogger::Instance()->SetRemark(false);
    JetScapeLogger::Instance()->SetVerboseLevel(0);
    Pythia8::Pythia pythia ("",false);//to be used for its particle info

    string directory = argv[1];
    //int RHIC_or_LHC = std::stoi(argv[2]);
    string input_filename = directory + "evolution_result.dat";

    vector<double> jet_radii = {0.2, 0.3, 0.4, 0.5, 0.6};
    size_t num_jet_radii = jet_radii.size();

    // Criteria of hadrons in Jet-Reconstruction.
    double hadron_eta_cut = 5;
    // if flag == 0, it's RHIC, otherwise it's LHC
    // and the pT_trk_min is 0.2 GeV if RHIC and 1 GeV otherwise
    //double pT_trk_min     = RHIC_or_LHC == 0 ? 0.2 : 1.0 ;//GeV

    auto reader=make_shared<JetScapeReaderAscii>(input_filename);

    // create the output files: particles and hadrons:
    std::ofstream fout (directory + "FinalStateParticles.txt");
    std::ofstream hout (directory + "FinalStateHadrons.txt");
    vector<std::ofstream> jet_files_ncut;
    //vector<std::ofstream> jet_files_wcut;
    //vector<std::ofstream> jet_files_chgd;
    jet_files_ncut.reserve(num_jet_radii);
    //jet_files_wcut.reserve(num_jet_radii);
    //jet_files_chgd.reserve(num_jet_radii);

    for (auto r: jet_radii)
    {
        string fname1 = directory + "FinalStateJets_" + to_string(r) + "_ncut.dat";
        //string fname2 = directory + "FinalStateJets_" + to_string(r) + "_wcut.dat";
        //string fname3 = directory + "FinalStateJets_" + to_string(r) + "_chg.dat";

        jet_files_ncut.emplace_back(ofstream(fname1));
        //jet_files_wcut.emplace_back(ofstream(fname2));
        //jet_files_chgd.emplace_back(ofstream(fname3));
    }

    vector<shared_ptr<Hadron>> hadrons;

    while (!reader->Finished())
    {
        reader->Next();
        fout <<"#Event_ID= "<< reader->GetCurrentEvent()+1 <<endl;
        hout <<"#Event_ID= "<< reader->GetCurrentEvent()+1 <<endl;
        for (int i_r = 0; i_r < num_jet_radii; i_r++)
        {
            jet_files_ncut[i_r] << "#Event_ID= " << reader->GetCurrentEvent()+1 << endl;            
            //jet_files_wcut[i_r] << "#Event_ID= " << reader->GetCurrentEvent()+1 << endl;
            //jet_files_chgd[i_r] << "#Event_ID= " << reader->GetCurrentEvent()+1 << endl;
        }

        // Write Showers to file:
        auto mShowers=reader->GetPartonShowers();
        fout <<"#Shower_Number= "<<mShowers.size()<<endl;
        for (int i=0;i<mShowers.size();i++)
        {
            fout<<"#Shower_ID= "<< i <<"\tParton_Number= "<<mShowers[i]->GetFinalPartons().size()<<endl;
            double shower_en = 0.0, shower_px = 0.0, shower_py = 0.0, shower_pz = 0.0;
            for ( int ipart = 0; ipart< mShowers[i]->GetFinalPartons().size(); ++ipart)
            {
                Parton p = *mShowers[i]->GetFinalPartons().at(ipart);
                fout <<"P\t"<< i <<","<< ipart <<","<< p.pid() <<","<< p.pstat() 
                     <<","  << p.e() << "," << p.pt() << "," << p.phi()<< ","<< p.eta() << endl;
                shower_en += p.e();
                shower_px += p.px();
                shower_py += p.py();
                shower_pz += p.pz();
            }
            double shower_pT  = sqrt(shower_px*shower_px + shower_py*shower_py);
            double shower_phi = atan2(shower_py, shower_px);
            double shower_eta = (fabs(shower_pz)>=fabs(shower_en)) ? 100. : atanh(shower_pz/shower_en);
            fout<<"S\t"<< i <<","<<shower_en<<","<<shower_pT<<","<<shower_phi<<","<<shower_eta<<endl;
        }

        // for this event, read in the hadrons
        // and write them to file.
        hadrons = reader->GetHadrons();
        fout <<"#Hadron_Number= "<< hadrons.size()<<endl;
        // populate these vectors here. they're independent of the jet cone size
        // so iterate over the hadron list once, populate these and then use the jet radius 
        // loop to construct jets
        vector<fjcore::PseudoJet> forFJ_ncut;
        //vector<fjcore::PseudoJet> forFJ_wcut;
        //vector<fjcore::PseudoJet> forFJ_chgd;
        //vector<fjcore::PseudoJet> antihadron;// to hold the holes, for later subtraction
        //bool charged;
        for (unsigned int i=0; i<hadrons.size(); i++)
        {
	        double hpsj_px = hadrons[i].get()->px(), hpsj_py = hadrons[i].get()->py();
	        double hpsj_pz = hadrons[i].get()->pz(), hpsj_e = hadrons[i].get()->e();

	        hout <<i<<","<< hadrons[i].get()->pstat()<<","<< hadrons[i].get()->pid()
	                <<","<< hadrons[i].get()->e()    <<","<< hadrons[i].get()->pt()
	                <<","<< hadrons[i].get()->phi()  <<","<< hadrons[i].get()->eta()<<endl;
            fjcore::PseudoJet Hadron_as_PseudoJet(hpsj_px, hpsj_py, hpsj_pz, hpsj_e);
	        Hadron_as_PseudoJet.set_user_index(hadrons[i].get()->pid());//User_index = pid

            //charged = fabs(pythia.particleData.charge(hadrons[i].get()->pid())) < 1e-2 ? false : true;
	        if ((hadrons[i].get()->pstat())==0)
	        {
	            //first if: just check if we're in the correct eta window
	            if ( fabs(hadrons[i].get()->eta()) < hadron_eta_cut)
	            {
	                forFJ_ncut.push_back(Hadron_as_PseudoJet);

	                //if (charged) 
                    //{
                    //    forFJ_chgd.push_back(Hadron_as_PseudoJet);
                    //    if ((hadrons[i].get()->pt())>pT_trk_min)
	                //        forFJ_wcut.push_back(Hadron_as_PseudoJet);
                    //}
                    //if ( not charged)
                    //    forFJ_wcut.push_back(Hadron_as_PseudoJet);
	            }
	        }
	        //else
	        //    antihadron.push_back(Hadron_as_PseudoJet);

        }
        
        //Now handle jets by looping over the the jet radii
        for (auto r : jet_radii)
        {
            fjcore::JetDefinition jet_def(fjcore::antikt_algorithm, r);// jet definition
            auto i_r = distance(jet_radii.begin(), find(jet_radii.begin(), jet_radii.end(), r));
            // cluster the jets
            fjcore::ClusterSequence seq_ncut(forFJ_ncut, jet_def);
            //fjcore::ClusterSequence seq_wcut(forFJ_wcut, jet_def);
            //fjcore::ClusterSequence seq_chgd(forFJ_chgd, jet_def);

            vector<fjcore::PseudoJet> jets_ncut = fjcore::sorted_by_pt(seq_ncut.inclusive_jets(10.));
            //vector<fjcore::PseudoJet> jets_wcut = fjcore::sorted_by_pt(seq_wcut.inclusive_jets(10.));
            //vector<fjcore::PseudoJet> jets_chgd = fjcore::sorted_by_pt(seq_chgd.inclusive_jets(10.));

            //loop over each of the above, if a hole falls inside a jet, subtract it.
            auto size_jets_ncut = jets_ncut.size();
            for (size_t i_jet = 0; i_jet < size_jets_ncut; i_jet++)
            {
                auto jet = jets_ncut[i_jet];
                double Jet_En  = jet.E();
                double Jet_Px  = jet.px();
                double Jet_Py  = jet.py();
                double Jet_Pz  = jet.pz();
                double Jet_eta = atanh(Jet_Pz/Jet_En);
                double Jet_pT  = sqrt(Jet_Px*Jet_Px + Jet_Py*Jet_Py);
                double Jet_phi = atan2(Jet_Py, Jet_Px);
                jet_files_ncut[i_r] << "J\t" << i_jet <<","<< Jet_En <<","<< Jet_pT <<","<< Jet_phi <<","<< Jet_eta <<endl;
                vector<fjcore::PseudoJet> hadrons_in_jet = fjcore::sorted_by_pt(jet.constituents());
                for (int ih = 0; ih < hadrons_in_jet.size(); ih++)
                {
                    jet_files_ncut[i_r] << "JH\t" << i_jet <<","<<ih<<","<<hadrons_in_jet[ih].user_index()
                      <<","<< hadrons_in_jet[ih].E()   <<"," << hadrons_in_jet[ih].pt()
                      <<","<< hadrons_in_jet[ih].phi() <<"," << hadrons_in_jet[ih].rap()<<endl;
                }
            }
            //auto size_jets_wcut = jets_wcut.size();
            //for (size_t i_jet = 0; i_jet < size_jets_wcut; i_jet++)
            //{
            //    auto jet = jets_wcut[i_jet];
            //    vector<PseudoJet> copy_antihadron(antihadron);
            //    auto it = copy_antihadron.begin();
            //    while (it != copy_antihadron.end())
            //    {
            //        // remove hole(negative) particle from jet
            //        double delta_r = it->delta_R(jet);
            //        if (delta_r < r)
            //        {
            //            jet -= *it;
            //            it = copy_antihadron.erase(it);
            //        }
            //        else it++;
            //    }
            //    //all possible holes removed from this jet. write it 
            //    //to the appropriate file
            //    double Jet_En  = jet.E();
            //    double Jet_Px  = jet.px();
            //    double Jet_Py  = jet.py();
            //    double Jet_Pz  = jet.pz();
            //    double Jet_eta = atanh(Jet_Pz/Jet_En);
            //    double Jet_pT  = sqrt(Jet_Px*Jet_Px + Jet_Py*Jet_Py);
            //    double Jet_phi = atan2(Jet_Py, Jet_Px);
            //    jet_files_wcut[i_r] << "J\t" << i_jet <<","<< Jet_En <<","<< Jet_pT <<","<< Jet_phi <<","<< Jet_eta <<endl;
            //    vector<fjcore::PseudoJet> hadrons_in_jet = fjcore::sorted_by_pt(jet.constituents());
            //    for (int ih = 0; ih < hadrons_in_jet.size(); ih++)
            //    {
            //        jet_files_wcut[i_r] << "JH\t" << i_jet <<","<<ih<<","<<hadrons_in_jet[ih].user_index()
            //          <<","<< hadrons_in_jet[ih].E()   <<"," << hadrons_in_jet[ih].pt()
            //          <<","<< hadrons_in_jet[ih].phi() <<"," << hadrons_in_jet[ih].rap()<<endl;
            //    }
            //}
            //auto size_jets_chgd = jets_chgd.size();
            //for (size_t i_jet = 0; i_jet < size_jets_chgd; i_jet++)
            //{
            //    auto jet = jets_chgd[i_jet];
            //    vector<PseudoJet> copy_antihadron(antihadron);
            //    auto it = copy_antihadron.begin();
            //    while (it != copy_antihadron.end())
            //    {
            //        // remove hole(negative) particle from jet
            //        double delta_r = it->delta_R(jet);
            //        if (delta_r < r)
            //        {
            //            jet -= *it;
            //            it = copy_antihadron.erase(it);
            //        }
            //        else it++;
            //    }
            //    //all possible holes removed from this jet. write it 
            //    //to the appropriate file
            //    double Jet_En  = jet.E();
            //    double Jet_Px  = jet.px();
            //    double Jet_Py  = jet.py();
            //    double Jet_Pz  = jet.pz();
            //    double Jet_eta = atanh(Jet_Pz/Jet_En);
            //    double Jet_pT  = sqrt(Jet_Px*Jet_Px + Jet_Py*Jet_Py);
            //    double Jet_phi = atan2(Jet_Py, Jet_Px);
            //    jet_files_chgd[i_r] << "J\t" << i_jet <<","<< Jet_En <<","<< Jet_pT <<","<< Jet_phi <<","<< Jet_eta <<endl;
            //    vector<fjcore::PseudoJet> hadrons_in_jet = fjcore::sorted_by_pt(jet.constituents());
            //    for (int ih = 0; ih < hadrons_in_jet.size(); ih++)
            //    {
            //        jet_files_chgd[i_r] << "JH\t" << i_jet <<","<<ih<<","<<hadrons_in_jet[ih].user_index()
            //          <<","<< hadrons_in_jet[ih].E()   <<"," << hadrons_in_jet[ih].pt()
            //          <<","<< hadrons_in_jet[ih].phi() <<"," << hadrons_in_jet[ih].rap()<<endl;
            //    }
            //}
        }
    }
    reader->Close();
    fout.close();
    hout.close();
    for (int i=0; i<num_jet_radii; i++)
    {
        jet_files_ncut[i].close();
        //jet_files_wcut[i].close();
        //jet_files_chgd[i].close();
    }
    cout << "Done with clustering and writing to file." << endl;
    return 0;
}
