#include "analysis.h"

void bin_partons_by_z(vector<Parton>* plist, vector<double>& gluon_hist, vector<double>& fermion_hist, vector<double>& zbins, double pTRef, double tau){
    size_t num_partons = plist->size();
    size_t num_zbins = zbins.size()-1;
    double z, weight;
    int id, index;
    for (size_t ih = 0; ih < num_partons; ih++){
        auto parton = plist->at(ih);
        if (abs(parton.p().eta()) > PARTON_ETA_CUT)
            continue;
        z = parton.p().pT()/pTRef;
        weight = parton.recoil() != -1.0 ? 1.0 : -1.0;
        id = abs(parton.id());
        if (parton.getShowerCreationTime() * 0.1973269804 > tau)
            continue;

        index = -1;
        for (size_t i = 0; i < num_zbins; i++){
            if ( zbins[i] <= z and z < zbins[i+1]){
                index = i;
                break;
            }
        }
        if (index < 0)
            continue;
        if (id <= 4)
            fermion_hist[index] += weight;
        else if (id == 21)
            gluon_hist[index] += weight;
        else
            continue;
    };

}
void bin_parton(vector<Parton>* plist, vector<double>& hist, vector<double>& pTbins, int recoil_flag, int particle_type){
    size_t num_partons = plist->size();
    size_t num_pTbins = pTbins.size()-1;
    double pT;
    int recoil;
    if (particle_type == QUARK){
        for (size_t ih = 0; ih < num_partons; ih++){
            auto parton = plist->at(ih);
            recoil = parton.recoil();
            if (abs(parton.id()) <= 4 and recoil == recoil_flag){
                if (fabs(parton.p().eta()) < CH_HAD_ETA_CUT){
                    pT = parton.p().pT();
                    for (size_t i=0; i < num_pTbins; i++){
                        if ( pTbins[i] <= pT and pT < pTbins[i+1]){
                            if (recoil != -1)
                                hist[i] += 1;
                            else
                                hist[i] -= 1;
                            break;
                        }
                    }
                }
            }
        }
    }
    else if (particle_type == GLUON){
        for (size_t ih = 0; ih < num_partons; ih++){
            auto parton = plist->at(ih);
            recoil = parton.recoil();
            if (parton.id() ==  21 and recoil == recoil_flag){
                if (fabs(parton.p().eta()) < CH_HAD_ETA_CUT){
                    pT = parton.p().pT();
                    for (size_t i=0; i < num_pTbins; i++){
                        if ( pTbins[i] <= pT and pT < pTbins[i+1]){
                            if (recoil != HOLE)
                                hist[i] += 1;
                            else
                                hist[i] -= 1;
                            break;
                        }
                    }
                }
            }
        }
    }
    else
        cout << "Particle Type flag not recognized: " << particle_type << endl;
}
void bin_chgd_hads(Event& evt, vector<double>& hist, vector<double>& pTbins, double weight){
    size_t num_hads_in_evt = evt.size();
    size_t num_pTbins = pTbins.size()-1;
    double pT;
    for (size_t ih = 0; ih < num_hads_in_evt; ih++){
        auto had = evt[ih];
        if (had.isHadron() and had.isFinal() and had.isCharged()){
            if (fabs(had.eta()) < CH_HAD_ETA_CUT){
                pT = had.pT();
                for (size_t i=0; i < num_pTbins; i++){
                    if ( pTbins[i] <= pT and pT < pTbins[i+1]){
                        hist[i] += weight;
                        break;
                    }
                }
            }
        }
    }
}

void bin_identified_had(Event& evt, vector<double>& hist, vector<double>& pTbins, double weight, int pid){
    size_t num_hads_in_evt = evt.size();
    size_t num_pTbins = pTbins.size()-1;
    double pT;
    int currID; 
    for (size_t ih = 0; ih < num_hads_in_evt; ih++){
        auto had = evt[ih];
        currID = abs(had.id());
        if (had.isHadron() and had.isFinal() and currID == pid){
            if (fabs(had.eta()) < CH_HAD_ETA_CUT){
                pT = had.pT();
                for (size_t i=0; i < num_pTbins; i++){
                    if ( pTbins[i] <= pT and pT < pTbins[i+1]){
                        hist[i] += weight;
                        break;
                    }
                }
            }
        }
    }
}
/*
 * Clustering Routines
 * */
void cluster_prep_partonic(vector<Parton>* plist, vector<PseudoJet>& ps, int hole_flag){
    size_t num_partons = plist->size();
    int id = 0;
    for (size_t ip = 0; ip < num_partons; ip++){
        auto part = plist->at(ip);
        id = part.id();
        if (id == 22)
            continue;
        if (fabs(part.p().eta()) < ETA_TRK_MAX and part.recoil() == hole_flag){
            double px = part.p().px();
            double py = part.p().py();
            double pz = part.p().pz();
            double ep = part.p().e();
            PseudoJet Hadron_as_PseudoJet(px, py, pz, ep);
            Hadron_as_PseudoJet.set_user_index(id);
            ps.push_back(Hadron_as_PseudoJet);
        }
    }
}

void cluster_prep_inclusive(Event& evt, vector<PseudoJet>& hads){
    size_t num_particles_in_evt = evt.size();
    for (size_t ip = 0; ip < num_particles_in_evt; ip++){
        auto part = evt[ip];
        if (part.isFinal() and fabs(part.p().eta()) < ETA_TRK_MAX){
            double px = part.p().px();
            double py = part.p().py();
            double pz = part.p().pz();
            double ep = part.p().e();
            PseudoJet Hadron_as_PseudoJet(px, py, pz, ep);
            if (fabs(part.charge()) > 0)
                Hadron_as_PseudoJet.set_user_index(1);
            else
                Hadron_as_PseudoJet.set_user_index(0);
            hads.push_back(Hadron_as_PseudoJet);
        }
    }
}

void cluster_prep_charged(Event& evt, vector<PseudoJet>& hads){
    size_t num_particles_in_evt = evt.size();
    for (size_t ip = 0; ip < num_particles_in_evt; ip++){
        auto part = evt[ip];
        if (part.isFinal() and fabs(part.eta()) < ETA_TRK_MAX and part.isHadron()){
            double px = part.px();
            double py = part.py();
            double pz = part.pz();
            double ep = part.eCalc();
            PseudoJet Hadron_as_PseudoJet(px, py, pz, ep);
            if (fabs(part.charge()) > 0)
                Hadron_as_PseudoJet.set_user_index(1);
            else
                Hadron_as_PseudoJet.set_user_index(0);
            hads.push_back(Hadron_as_PseudoJet);
        }
    }
}

void cluster_prep_inclusive_wth_trk_cut(Event& evt, vector<PseudoJet>& hads, double trk_cut){
    size_t num_particles_in_evt = evt.size();
    for (size_t ip = 0; ip < num_particles_in_evt; ip++){
        auto part = evt[ip];
        if (part.isFinal() and fabs(part.eta()) < ETA_TRK_MAX){
            double px, py, pz, ep;
            if (part.isCharged() and part.pT() > trk_cut){
                px = part.px();
                py = part.py();
                pz = part.pz();
                ep = part.eCalc();
            }
            else{
                px = part.px();
                py = part.py();
                pz = part.pz();
                ep = part.eCalc();
            }
            PseudoJet Hadron_as_PseudoJet(px, py, pz, ep);
            if (fabs(part.charge()) > 0)
                Hadron_as_PseudoJet.set_user_index(1);
            else
                Hadron_as_PseudoJet.set_user_index(0);
            hads.push_back(Hadron_as_PseudoJet);
        }
    }
}

/*
 * Function to subtract holes from jets
 * */
void subtract_holes_from_jets(vector<PseudoJet>& jets, vector<PseudoJet>& holes, double R){
    auto num_jets = jets.size();
    double delta_r;
    for (size_t ij = 0; ij < num_jets; ij++){
        auto jet = jets[ij];
        vector<PseudoJet> copy_holes(holes);
        auto it = copy_holes.begin();
        while (it != copy_holes.end()){
            // remove hole(negative) particle from jet
            delta_r = it->delta_R(jet);
            if (delta_r < R){
                jet -= *it;
                it = copy_holes.erase(it);
            }
            else it++;
        }
    }
}

void bin_jets(vector<PseudoJet>& jets, vector<double>& hist, vector<double>& pTbins){
    size_t num_jets = jets.size();
    size_t num_pTbins = pTbins.size()-1;
    double pT;
    for (size_t ij = 0; ij < num_jets; ij++){
        auto jet = jets[ij];
        // calculate jet eta:
        if (fabs(jet.pseudorapidity()) < JET_ETA_CUT){
            pT = jet.perp();
            for (size_t i = 0; i < num_pTbins; i++){
                if (pTbins[i] <= pT and pT < pTbins[i+1]){
                    hist[i] += 1;
                    break;
                }
            }
        }
    }
}

int bin_jet_fragmentation(vector<PseudoJet>& jets, vector<double>& FFz, 
     vector<double>& zbins, vector<double>& FFpT, vector<double>& pTbins){
   
    int num_jets_processed = 0;
    size_t num_jets = jets.size();
    size_t num_zbins = zbins.size()-1;
    size_t num_ptbins = pTbins.size()-1;

    double jet_eta, jet_pT, jet_pT_squared;
    double jet_en, hadron_en;
    double jet_phi, hadron_phi;
    double hadron_pT, hadron_z, hadron_eta;
    int hadron_charge;
    double delta_phi;
    double pJet_dot_pJet, pJet_dot_pTrk;
    double z;
    
    for (size_t ijet = 0; ijet < num_jets; ijet++){
        auto jet = jets[ijet];

        jet_pT = jet.pt();
        jet_eta = jet.pseudorapidity();

        if (jet_pT < JET_PTMIN_FF  or jet_pT > JET_PTMAX_FF)
            continue;
        if (fabs(jet_eta) > JET_ETA_FF)
            continue;
        num_jets_processed += 1;
        jet_en = jet.e();
        jet_pT_squared = jet_pT*jet_pT;
        pJet_dot_pJet = jet_pT_squared + pow(jet_en*tanh(jet_eta), 2);
        for (auto hadron : jet.constituents()){
            hadron_charge = hadron.user_index();

            if ( not hadron_charge ) continue;

            hadron_pT = hadron.pt();
            hadron_eta = hadron.pseudorapidity();
            hadron_en = hadron.e();


            delta_phi = jet.delta_phi_to(hadron);

            pJet_dot_pTrk = jet_pT*hadron_pT*cos(delta_phi) + jet_en*hadron_en*tanh(jet_eta)*tanh(hadron_eta);
            z = pJet_dot_pTrk / pJet_dot_pJet;

            for (size_t iz = 0; iz < num_zbins; iz++){
                if (zbins[iz] <= z and z <zbins[iz+1]){
                    FFz[iz] += 1;
                    break;
                }
            }
            // Now bin the pT
            for (size_t ipt=0; ipt < num_ptbins; ipt++){
                if ( pTbins[ipt]<= hadron_pT and hadron_pT < pTbins[ipt+1]){
                    FFpT[ipt] += 1;
                    break;
                }
            }
        }
    }
    return num_jets_processed;
}

int bin_jet_shape(vector<PseudoJet>& jets, vector<double>& shape, vector<double>& rbins){
    int num_processed_jets = 0;
    size_t num_jets = jets.size();
    size_t num_rbins = rbins.size()-1;
    double jet_pt, jet_eta;
    double hadron_pt, hadron_eta;
    int hadron_charge;
    double delta_phi, delta_eta, r;
    for (size_t ij=0; ij < num_jets; ij++){
        auto jet = jets[ij];
        jet_pt = jet.pt();
        jet_eta = jet.pseudorapidity();
        if (jet_pt < JET_PTMIN_SHAPE or fabs(jet_eta) < JET_ETA_MIN_SHAPE or fabs(jet_eta) > JET_ETA_MAX_SHAPE)
            continue;
        num_processed_jets += 1;
        for (auto hadron : jet.constituents()){
            hadron_charge = hadron.user_index();
            if (not hadron_charge) continue;
            hadron_pt = hadron.pt();
            hadron_eta = hadron.pseudorapidity();
            delta_phi = jet.delta_phi_to(hadron);
            delta_eta = jet_eta - hadron_eta;
            r = sqrt(delta_phi*delta_phi + delta_eta*delta_eta);
            // now bin it with weight hadron_pt/jet_pt
            for (size_t ir=0; ir < num_rbins; ir++){
                if (rbins[ir] <= r and r < rbins[ir+1]){
                    shape[ir] += (hadron_pt/jet_pt);
                    break;
                }
            }
        }
    }
    return num_processed_jets;
}

int bin_partonic_shape(vector<PseudoJet>& jets, vector<double>& shape, vector<double>& rbins, int particle_type){
    int num_processed_jets = 0;
    size_t num_jets = jets.size();
    size_t num_rbins = rbins.size()-1;
    double jet_pt, jet_eta;
    double parton_pt, parton_eta;
    double delta_phi, delta_eta, r;
    int parton_id, my_type;
    for (size_t ij=0; ij < num_jets; ij++){
        auto jet = jets[ij];
        jet_pt = jet.pt();
        jet_eta = jet.pseudorapidity();
        if (jet_pt < PARTON_JET_PTMIN_SHAPE or fabs(jet_eta) < PARTON_JET_ETA_MIN_SHAPE or fabs(jet_eta) > PARTON_JET_ETA_MAX_SHAPE)
            continue;
        num_processed_jets += 1;
        for (auto parton : jet.constituents()){
            parton_id = abs(parton.user_index());
            my_type = parton_id < 4 ? QUARK : GLUON;
            if (my_type != particle_type ) continue;
            parton_pt  = parton.pt();
            parton_eta = parton.pseudorapidity();
            delta_phi = jet.delta_phi_to(parton);
            delta_eta = jet_eta - parton_eta;
            r = sqrt(delta_phi*delta_phi + delta_eta*delta_eta);
            // now bin it with weight parton_pt/jet_pt
            for (size_t ir=0; ir < num_rbins; ir++){
                if (rbins[ir] <= r and r < rbins[ir+1]){
                    shape[ir] += (parton_pt/jet_pt);
                    break;
                }
            }
        }
    }
    return num_processed_jets;
}
void bin_photons(vector<Parton>* plist, vector<double>& prompt, vector<double>& conv, vector<double>& brem, vector<double>& pTbins){
    size_t num_partons = plist->size();
    size_t num_pTbins = pTbins.size()-1;
    double pt,eta;
    int source, id;

    for (size_t ip = 0; ip < num_partons; ip++){
        auto parton = plist->at(ip);
        id = parton.id();
        if (id != 22)//not a photon
            continue;
        eta = fabs(parton.p().eta());
        if (eta > PHOTON_ETA_CUT)
            continue;
        pt = parton.p().pT();
        source = parton.source();
        for (size_t ibin=0; ibin < num_pTbins; ibin++){
            if (pTbins[ibin] <= pt and pt < pTbins[ibin+1]){
                if ( source == 0)
                    prompt[ibin] += 1;
                else if (source == 1)
                    conv[ibin] += 1;
                else if (source == 2)
                    brem[ibin] += 1;
                else
                    cout << "Found a photon with a strange source: " << source << endl;
                break;
            }
        }
    }
}
