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

    string input_filename = "cujet_single_out.dat";
    // Criteria of hadrons in Jet-Reconstruction.
    double Jet_Radius = 0.4;
    double pT_min = 0; //GeV
    double eta_min = -5;
    double eta_max = 5;
    if ( argc >1 && string(argv[1]) == "-h" ) {
        cout << "Usage: ./cujet_reader (input_filename) (pT_min) (eta_max) (eta_min)" << endl;
        return -1;
    }
    stringstream ss_eta_min, ss_eta_max, ss_pt;
    switch ( argc ){
        break;
    case 5:
        ss_eta_min << string(argv[4]);
        ss_eta_min >> eta_min;
    case 4:
        ss_eta_max << string(argv[3]);
        ss_eta_max >> eta_max;
    case 3:
        ss_pt << string(argv[2]);
        ss_pt >> pT_min;
    case 2:
        input_filename=string(argv[1]);
        break;
    case 1:
        break;
    case 0:
        break;
    default:
        cout << "Usage: cujet_reader [input_filename]" << endl;
        return -1;
    }
    if (argc<5) { eta_min = -eta_max; }
    auto reader=make_shared<JetScapeReaderAscii>(input_filename);
    std::ofstream fout ("FinalStateParticles.txt");
    std::ofstream hout ("FinalStateHadrons.txt");
    std::ofstream jout ("FinalStateJets.txt");
    vector<shared_ptr<Hadron>> hadrons;
    fjcore::JetDefinition jet_def(fjcore::antikt_algorithm, Jet_Radius);
    while (!reader->Finished())
    {
        reader->Next();
        fout <<"#Event_ID= "<< reader->GetCurrentEvent()+1 <<endl;
        hout <<"#Event_ID= "<< reader->GetCurrentEvent()+1 <<endl;
        jout <<"#Event_ID= "<< reader->GetCurrentEvent()+1 <<endl;

        auto mShowers=reader->GetPartonShowers();
        fout <<"#Shower_Number= "<<mShowers.size()<<endl;
        for (int i=0;i<mShowers.size();i++)
        {
            fout<<"#Shower_ID= "<< i <<"\tParton_Number= "<<mShowers[i]->GetFinalPartons().size()<<endl;
            double shower_en = 0.0, shower_px = 0.0, shower_py = 0.0, shower_pz = 0.0;
            for ( int ipart = 0; ipart< mShowers[i]->GetFinalPartons().size(); ++ipart)
            {
                Parton p = *mShowers[i]->GetFinalPartons().at(ipart);
                //if(abs(p.pid())!=5) continue;
                fout <<"P\t"<< i <<"\t"<< ipart <<"\t"<< p.pid() //<<"\t"<< p.pstat() 
                  << "\t"<< p.e() << "\t" << p.pt() << "\t" << p.phi()<< "\t"<< p.eta() << endl;
                shower_en += p.e();
                shower_px += p.px();
                shower_py += p.py();
                shower_pz += p.pz();
            }
            double shower_pT  = sqrt(shower_px*shower_px + shower_py*shower_py);
            double shower_phi = atan2(shower_py, shower_px);
            double shower_eta = (fabs(shower_pz)>=fabs(shower_en)) ? 100. : atanh(shower_pz/shower_en);
            fout<<"S\t"<< i <<"\t"<<shower_en<<"\t"<<shower_pT<<"\t"<<shower_phi<<"\t"<<shower_eta<<endl;
        }

        vector<fjcore::PseudoJet> forFJ;
        vector<fjcore::PseudoJet> sbtrHadron;
        hadrons = reader->GetHadrons();
        fout <<"#Hadron_Number= "<< hadrons.size()<<endl;
        for(unsigned int i=0; i<hadrons.size(); i++)
        {
            double hpsj_px = hadrons[i].get()->px(), hpsj_py = hadrons[i].get()->py();
            double hpsj_pz = hadrons[i].get()->pz(), hpsj_e = hadrons[i].get()->e();
            fjcore::PseudoJet Hadron_as_PseudoJet(hpsj_px, hpsj_py, hpsj_pz, hpsj_e);
            Hadron_as_PseudoJet.set_user_index(hadrons[i].get()->pid());//User_index = pid
            if ((hadrons[i].get()->pstat())==0)
            {
                if( ((hadrons[i].get()->pt())>pT_min)
                    &&((hadrons[i].get()->eta())>eta_min)
                    &&((hadrons[i].get()->eta())<eta_max)
                   )
                {
                    forFJ.push_back(Hadron_as_PseudoJet);
                }
            }
            else
            {
                sbtrHadron.push_back(Hadron_as_PseudoJet);
            }
            fout<<"H\t"<<i<<"\t"<<hadrons[i].get()->pstat()<<"\t"<<hadrons[i].get()->pid()
              <<"\t"<< hadrons[i].get()->e() <<"\t"<< hadrons[i].get()->pt()
              <<"\t"<< hadrons[i].get()->phi() <<"\t"<< hadrons[i].get()->eta()<<endl;
            hout<<i<<"\t"<<hadrons[i].get()->pstat()<<"\t"<<hadrons[i].get()->pid()
              <<"\t"<< hadrons[i].get()->e() <<"\t"<< hadrons[i].get()->pt()
              <<"\t"<< hadrons[i].get()->phi() <<"\t"<< hadrons[i].get()->eta()<<endl;
        }
        fjcore::ClusterSequence hcs(forFJ, jet_def);
        vector<fjcore::PseudoJet> hjets = fjcore::sorted_by_pt(hcs.inclusive_jets());
        fout<<"#Anti-kT_Jet_Number= " << hjets.size() << endl;
        for (int k=0;k<hjets.size();k++)
        {
            double Jet_En = hjets[k].E();
            double Jet_Px = hjets[k].px();
            double Jet_Py = hjets[k].py();
            double Jet_Pz = hjets[k].pz();
            for (int ah=0; ah<sbtrHadron.size(); ah++)
            {
                if ((sbtrHadron[ah].delta_R(hjets[k]) < Jet_Radius)
                    &&(sbtrHadron[ah].eta()<eta_max)
                    &&(sbtrHadron[ah].eta()>eta_min))
                {
                    Jet_En -= sbtrHadron[ah].E();
                    Jet_Px -= sbtrHadron[ah].px();
                    Jet_Py -= sbtrHadron[ah].py();
                    Jet_Pz -= sbtrHadron[ah].pz();
                }
            }
            double Jet_eta = atanh(Jet_Pz/Jet_En);
            double Jet_pT  = sqrt(Jet_Px*Jet_Px + Jet_Py*Jet_Py);
            double Jet_phi = atan2(Jet_Py, Jet_Px);
            jout << k <<"\t"<< Jet_En <<"\t"<< Jet_pT <<"\t"<< Jet_phi <<"\t"<< Jet_eta <<endl;
            fout <<"J\t"<< k <<"\t"<<hjets[k].constituents().size()<<"\t"<< Jet_En
              <<"\t"<< Jet_pT <<"\t"<< Jet_phi <<"\t"<< Jet_eta <<endl;
            fout <<"Junc\t"<< k <<"\t"<<hjets[k].constituents().size()<<"\t"<< hjets[k].E() 
              <<"\t" << hjets[k].pt() <<"\t"<< hjets[k].phi() <<"\t"<< hjets[k].rap() <<endl;
            
            vector<fjcore::PseudoJet> hadron_in_jets = fjcore::sorted_by_pt(hjets[k].constituents());
            for(int l=0; l<hadron_in_jets.size(); l++)
            {
                fout<<"JH\t"<<k<<"\t"<<l<<"\t"<<hadron_in_jets[l].user_index()
                  <<"\t"<< hadron_in_jets[l].E() <<"\t" << hadron_in_jets[l].pt()
                  <<"\t"<< hadron_in_jets[l].phi() <<"\t"<< hadron_in_jets[l].rap()<<endl;
            }
        }

    }
    reader->Close();
    fout.close();
    hout.close();
    jout.close();
}
