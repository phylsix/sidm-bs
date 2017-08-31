// -*- C++ -*-
//
// Package:    sidm-bs/mcValidation
// Class:      mcValidation
#include <algorithm>
#include <cmath>

#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/Common/interface/Ptr.h"
#include "DataFormats/Math/interface/LorentzVector.h"
#include "DataFormats/Math/interface/deltaR.h"

#include "sidm-bs/jetComponent/interface/physicsObject.h"
#include "sidm-bs/jetComponent/interface/mcValidation.h"
#include "sidm-bs/jetComponent/interface/utilities.h"

sidm::mcValidation::mcValidation(const edm::ParameterSet& iConfig):
    genParticleTk_(consumes<edm::View<reco::GenParticle> >(iConfig.getUntrackedParameter<edm::InputTag>("GenParticleTag_", edm::InputTag("prunedGenParticles")))),
    ssVerticeTk_(consumes<edm::View<reco::VertexCompositePtrCandidate> >(iConfig.getUntrackedParameter<edm::InputTag>("SsVerticeTag_", edm::InputTag("slimmedSecondaryVertices")))),
    patElectronTk_(consumes<edm::View<pat::Electron> >(iConfig.getUntrackedParameter<edm::InputTag>("PatElectronTag_", edm::InputTag("slimmedElectrons")))),
    patJetTk_(consumes<edm::View<pat::Jet> >(iConfig.getUntrackedParameter<edm::InputTag>("PatJetTag_", edm::InputTag("slimmedJets")))),
    zpMassSb_(iConfig.getUntrackedParameter<double>("ZpMassSideBand")),
    zpMass_(iConfig.getUntrackedParameter<double>("ZpMass")),
    dRusb_(iConfig.getUntrackedParameter<double>("dRusb"))
{
    usesResource("TFileService");
    eventNum_ = 0;
    event2pairs_ = 0;
    event2pairsEpEp_ = 0;
    event2pairsEpEj_ = 0;
    event2pairsEjEj_ = 0;
    event2pairsEpJj_ = 0;
    event2pairsEjJj_ = 0;
    event2pairsJjJj_ = 0;
    event1pair_ = 0;
    event1pairEp_ = 0;
    event1pairEj_ = 0;
    event1pairJj_ = 0;
    event0pair_ = 0;
}


sidm::mcValidation::~mcValidation()
{

    std::cout << "\n\nNumber of events: " << eventNum_ << "\n\n";
    std::cout << std::setprecision(4) << std::fixed;
    std::cout << "2 pairs events: " << event2pairs_ << " [" << float(event2pairs_)/eventNum_ << "]\n";
    std::cout << "\tEpEp: " << std::setw(7) << event2pairsEpEp_ << " [" << float(event2pairsEpEp_)/event2pairs_ << "][" << float(event2pairsEpEp_)/eventNum_ <<"]\n";
    std::cout << "\tEpEj: " << std::setw(7) << event2pairsEpEj_ << " [" << float(event2pairsEpEj_)/event2pairs_ << "][" << float(event2pairsEpEj_)/eventNum_ <<"]\n";
    std::cout << "\tEjEj: " << std::setw(7) << event2pairsEjEj_ << " [" << float(event2pairsEjEj_)/event2pairs_ << "][" << float(event2pairsEjEj_)/eventNum_ <<"]\n";
    std::cout << "\tEpJj: " << std::setw(7) << event2pairsEpJj_ << " [" << float(event2pairsEpJj_)/event2pairs_ << "][" << float(event2pairsEpJj_)/eventNum_ <<"]\n";
    std::cout << "\tEjJj: " << std::setw(7) << event2pairsEjJj_ << " [" << float(event2pairsEjJj_)/event2pairs_ << "][" << float(event2pairsEjJj_)/eventNum_ <<"]\n";
    std::cout << "\tJjJj: " << std::setw(7) << event2pairsJjJj_ << " [" << float(event2pairsJjJj_)/event2pairs_ << "][" << float(event2pairsJjJj_)/eventNum_ <<"]\n";
    std::cout << "1 pair events: " << event1pair_ << " [" << float(event1pair_)/eventNum_ << "]\n";
    std::cout << "\tEp: " << std::setw(7) << event1pairEp_ << " [" << float(event1pairEp_)/event1pair_ << "][" << float(event1pairEp_)/eventNum_ <<"]\n";
    std::cout << "\tEj: " << std::setw(7) << event1pairEj_ << " [" << float(event1pairEj_)/event1pair_ << "][" << float(event1pairEj_)/eventNum_ <<"]\n";
    std::cout << "\tEp: " << std::setw(7) << event1pairJj_ << " [" << float(event1pairJj_)/event1pair_ << "][" << float(event1pairJj_)/eventNum_ <<"]\n";
    std::cout << "0 pair events: " << event0pair_ << " [" << float(event0pair_)/eventNum_ << "]\n";
    eventNum_ = 0;

}



// ------------ method called for each event  ------------
void
sidm::mcValidation::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
    using namespace edm;
    using namespace std;

    Handle<View<reco::GenParticle> > genParticleHdl_;
    iEvent.getByToken(genParticleTk_, genParticleHdl_);

    Handle<View<reco::VertexCompositePtrCandidate> > ssVerticeHdl_;
    iEvent.getByToken(ssVerticeTk_, ssVerticeHdl_);

    Handle<View<pat::Electron> > patElectronHdl_;
    iEvent.getByToken(patElectronTk_, patElectronHdl_);

    Handle<View<pat::Jet> > patJetHdl_;
    iEvent.getByToken(patJetTk_, patJetHdl_);


    //--------------------
    //----GEN PARTICLE----
    //--------------------

    // Number of electrons in genParticle collection: 
    // 1. final states electron. 2. pt > 2GeV.
    electron_N = std::count_if(cbegin(*genParticleHdl_), cend(*genParticleHdl_),
                               [](const auto& p){ return p.pdgId() == 11 && 
                                                         p.status() == 1 &&
                                                         p.pt() > 2. ; });

    // Number of electrons whose mother is dark photon in genParticle collection.
    electron_from_zp_N = std::count_if(cbegin(*genParticleHdl_),cend(*genParticleHdl_),
                                       [](const auto& p){ return p.pdgId() == 11 &&
                                                                 p.mother()->pdgId() == 32; });

    // Number of positrons in genParticle collection:
    // 1. final states positron. 2. pt > 2GeV.
    positron_N = std::count_if(cbegin(*genParticleHdl_), cend(*genParticleHdl_),
                               [](const auto& p){ return p.pdgId() == -11 &&
                                                         p.status() == 1  &&
                                                         p.pt() > 2. ; });

    // Number of positrons whose mother is dark photon in genParticle collection.
    positron_from_zp_N = std::count_if(cbegin(*genParticleHdl_), cend(*genParticleHdl_),
                                       [](const auto& p){ return p.pdgId() == -11 &&
                                                                 p.mother()->pdgId() == 32; });

    // Number of dark photons in genParticle collection.
    zp_N = std::count_if(cbegin(*genParticleHdl_), cend(*genParticleHdl_),
                         [](const auto& p){ return p.pdgId() == 32; });

    // Number of pseudoscalar in genParticle collection.
    ps_N = std::count_if(cbegin(*genParticleHdl_), cend(*genParticleHdl_),
                         [](const auto& p){ return p.pdgId() == 35; });

    // Number of electrons in pat::Electron collection.
    patE_N = patElectronHdl_->size();

    // Number of jets in pat::Jet collection.
    patJet_N = patJetHdl_->size();

    eventTree_->Fill();

    //------------------------------------------
    // Store dark photon object for future matching.
    vector<sidm::Zp> genZps{};

    for ( const auto& p : *genParticleHdl_ ) {
        if ( abs(p.pdgId()) != 32 && 
             abs(p.pdgId()) != 35 ) {continue;} // Zp or Ps

        /// Zp
        if ( p.pdgId() == 32 ) {

            if ( p.numberOfDaughters() != 2 ) {
                cout << "[" << eventNum_ << "] "
                     << "Found a darkphoton OR pscalor whose daughter is NOT 2 !!" << endl;
                cout << "--> Instead: pdgId"<< p.pdgId() << " has " << p.numberOfDaughters() << endl;
                continue;
            }

            if ( p.daughter(0)->pdgId() * p.daughter(1)->pdgId() != -11*11 ) {
                cout << "[" << eventNum_ << "] "
                     << "Found a darkphoton whose daughters are not electron pair" << endl;
                cout << "--> Instead (" << p.daughter(0)->pdgId() << ", "
                                        << p.daughter(1)->pdgId() << ")" <<endl;
                continue;
            }

            if ( p.daughter(0)->pdgId() == 11 ) {
                // 0: electron, 1: positron
                zp_ = sidm::Zp(p.daughter(0), p.daughter(1));
                genZps.push_back(zp_);

            } else {
                // 0: positron, 1: electron
                zp_ = sidm::Zp(p.daughter(1), p.daughter(0));
                genZps.push_back(zp_);
            }

            zp_._eventId = eventNum_;
            zp_._pt      = p.pt();
            zp_._eta     = p.eta();
            zp_._mass    = p.mass();

            darkPhoton_reco_->Fill(); 

        }

        /// Ps
        if ( p.pdgId() == 35 ) {
        
            if ( p.numberOfDaughters() != 2 ) { 
                continue;
            }

            if ( p.daughter(0)->pdgId() * p.daughter(1)->pdgId() != 32*32 ) {
                cout << "[" << eventNum_ << "] "
                     << "Found a pscalar whose daughters are not dark photons," <<endl;
                cout << "--> Instead (" << p.daughter(0)->pdgId() << ", "
                                        << p.daughter(1)->pdgId() << ")" <<endl;
                continue;
            }

            ps_._eventId = eventNum_;
            ps_._mass    = p.mass();
            ps_._invM    = ( p.daughter(0)->p4() + p.daughter(1)->p4() ).M();
            ps_._dEta    = std::abs( p.daughter(0)->eta() - p.daughter(1)->eta() );
            ps_._dPhi    = std::abs( p.daughter(0)->phi() - p.daughter(1)->phi() );
            ps_._dR      = deltaR( *(p.daughter(0)), *(p.daughter(1)) );

            pscalar_reco_->Fill();

        }
    }

    //--------------------
    //----slimmedSecondaryVertices----
    //--------------------

    /*--
    cout << "[" << eventNum_ << "] "
         << "Num.Of.SecondaryVertices: " << ssVerticeHdl_->size() << endl;
    --*/


    //--------------------
    //----GEN MATCHING----
    //--------------------

    
    int PAIR(-1);                            // Number of reconstructed pairs.
    int EPPAIR(-1), EJPAIR(-1), JJPAIR(-1);  // Categorize them

    
    vector<Ptr<pat::Electron> > patElectronPtr_ = patElectronHdl_->ptrs();  // edm::Ptr can be used to compair and find.
    vector<Ptr<pat::Electron> > patEp_es{}, patEp_ps{};                     // store electron ptrs & positron ptrs separately.

    copy_if(patElectronPtr_.begin(), patElectronPtr_.end(), back_inserter(patEp_es),
            [](const Ptr<pat::Electron>& p){return p->charge()<0;});
    copy_if(patElectronPtr_.begin(), patElectronPtr_.end(), back_inserter(patEp_ps),
            [](const Ptr<pat::Electron>& p){return p->charge()>0;});

    /**
     * __EP PAIR SEARCH__
     */
    if ( patEp_es.size()>0 && patEp_ps.size()>0 ) {
        // At least a positron and an electrons are reconstructed
        // construct pair vector
        sidm::pairvec<pat::Electron> epPair(patEp_es, patEp_ps);
        vector<pair<Ptr<pat::Electron>, Ptr<pat::Electron> > > epPairTmp = epPair.get();

        // narrow mass sideband
        epPairTmp.erase(remove_if(begin(epPairTmp), end(epPairTmp), [this](const auto& p){ 
                    float m = (p.first->p4()+p.second->p4()).M();
                    return m<(1-zpMassSb_)*zpMass_ || m>(1+zpMassSb_)*zpMass_;
                    }), epPairTmp.end());

        if ( epPairTmp.size() == 0 ) {
            // no pairs in mass band..
            PAIR   = 0;
            EPPAIR = 0;
        } else if (epPairTmp.size()==1) {
            // ONLY 1 pair in mass band 
            // --> looking at dR matched or not with one of the two dark photons.
            sidm::match_patPair_with_zps(epPairTmp[0], genZps, dRusb_); //< 1. pat::Electron pair 2. gen dark photon vector. 3. dR limit
            
            if (count_if(cbegin(genZps), cend(genZps), [](const sidm::Zp& p){return p.matched;})) {
                PAIR   = 1; EPPAIR = 1;
                zp_r_ = sidm::Zp(epPairTmp[0]);
                zp_r_._eventId = eventNum_;
                darkPhoton_rereco_->Fill();
                // Remove the matched pat pair from pat::Electron collection.
                sidm::remove_from_collection(&patElectronPtr_, epPairTmp[0].first);
                sidm::remove_from_collection(&patElectronPtr_, epPairTmp[0].second);
            } else {
                PAIR   = 0;
                EPPAIR = 0;
            }
        } else if (epPairTmp.size()>1) {
            sort(begin(epPairTmp), end(epPairTmp), [this](const auto& lhs, const auto& rhs){
                    float lhsM = (lhs.first->p4()+lhs.second->p4()).M();
                    float rhsM = (rhs.first->p4()+rhs.second->p4()).M();
                    return std::abs(lhsM-zpMass_) < std::abs(rhsM-zpMass_);
                    });
            sidm::pairvec<pat::Electron> epPair(epPairTmp, true);
            epPairTmp.clear();
            epPairTmp = epPair.get_zip();

            if (epPairTmp.size() == 1) {
                /// ONLY 1 unique pair, esscencially the same way as we deal with 1 pair
                sidm::match_patPair_with_zps(epPairTmp[0], genZps, dRusb_); //< 1. pat::Electron pair 2. gen dark photon vector. 3. dR limit

                if (count_if(cbegin(genZps), cend(genZps), [](const sidm::Zp& p){return p.matched;})) {
                    PAIR = 1; EPPAIR = 1;
                    zp_r_ = sidm::Zp(epPairTmp[0]);
                    zp_r_._eventId = eventNum_;
                    darkPhoton_rereco_->Fill();
                    // Remove the matched pat pair from pat::Electron collection.
                    sidm::remove_from_collection(&patElectronPtr_, epPairTmp[0].first);
                    sidm::remove_from_collection(&patElectronPtr_, epPairTmp[0].second);
                } else { PAIR = 0; EPPAIR = 0; }
            } else {
                /// More than 1 unique pair in mass band.
                /// Keep record which zp has already been matched, skip if true
                bool first_genZp_matched(false);
                bool second_genZp_matched(false);
                for (const auto& q : epPairTmp) {
                    sidm::match_patPair_with_zps(q, genZps, dRusb_); //< This guarentee only 1 or 0 of them would be matched
                    if (!first_genZp_matched && genZps[0].matched) {
                        first_genZp_matched = true;
                        zp_r_ = sidm::Zp(q);
                        zp_r_._eventId = eventNum_;
                        darkPhoton_rereco_->Fill();
                        // Remove the matched pat pair from pat::Electron collection.
                        sidm::remove_from_collection(&patElectronPtr_, q.first);
                        sidm::remove_from_collection(&patElectronPtr_, q.second);
                    }
                    if (!second_genZp_matched && genZps[1].matched) {
                        second_genZp_matched = true;
                        zp_r_ = sidm::Zp(q);
                        zp_r_._eventId = eventNum_;
                        darkPhoton_rereco_->Fill();
                        // Remove the matched pat pair from pat::Electron collection.
                        sidm::remove_from_collection(&patElectronPtr_, q.first);
                        sidm::remove_from_collection(&patElectronPtr_, q.second);
                    }
                    /// If both matched, skip the rest.
                    if (first_genZp_matched && second_genZp_matched) break;
                }

                int matchedNum = static_cast<int>(first_genZp_matched) + static_cast<int>(second_genZp_matched);
                if (matchedNum == 2) {
                    PAIR = 2;           //< 2 Pairs are reconstructed.
                    EPPAIR = 2;         //< Both are electron-positron pairs.
                    EJPAIR = JJPAIR = 0;
                    ++event2pairs_;
                    ++event2pairsEpEp_;
                }
                else if (matchedNum == 1) { PAIR = 1; EPPAIR = 1; }
                else if (matchedNum == 0) { PAIR = 0; EPPAIR = 0; }
                
            }
        }
            
    } else {
        /**
         * 1. Only electrons/positrons in pat::Electron collection
         * 2. No electrons/positrons in pat::Electron collection
         * In neither situation electron pair would be reconstructed.
         */
        PAIR = 0; EPPAIR = 0;
    }
    // ---------------EP SEARCH END---------------


    /**
     * __EJ PAIR SEARCH__
     */
    vector<Ptr<pat::Jet> > patJetPtr_ = patJetHdl_->ptrs();
    if (PAIR < 2) {
        sidm::pairvec<pat::Electron, pat::Jet> ejPair(patElectronPtr_, patJetPtr_);
        vector<pair<Ptr<pat::Electron>, Ptr<pat::Jet> > > ejPairTmp = ejPair.get();

        // narrow mass sideband
        ejPairTmp.erase(remove_if(begin(ejPairTmp), end(ejPairTmp), [this](const auto& p) {
                    float m = (p.first->p4()+p.second->p4()).M();
                    return m<(1-zpMassSb_)*zpMass_ || m>(1+zpMassSb_)*zpMass_;
                    }), ejPairTmp.end());
        if (ejPairTmp.size() == 0) {
            // No pairs in mass band..
            EJPAIR = 0; //< PAIR remains as same (1).
        } else {
            // First, sort by invM
            sort(begin(ejPairTmp), end(ejPairTmp), [this](const auto& lhs, const auto& rhs) {
                    float lhsM = (lhs.first->p4()+lhs.second->p4()).M();
                    float rhsM = (rhs.first->p4()+rhs.second->p4()).M();
                    return std::abs(lhsM-zpMass_) < std::abs(rhsM-zpMass_);
                    });
            sidm::pairvec<pat::Electron, pat::Jet> ejPair(ejPairTmp, true);
            ejPairTmp.clear();
            ejPairTmp = ejPair.get_zip();

            if (PAIR == 1) {
                // We only need to look at the first one, aka, closest to the mass,
                // compare dR with one of the two gen darkPhotons who is not matched yet.
                for (const auto& q : ejPairTmp) {
                    sidm::match_patPair_with_zps(q, genZps, dRusb_);
                    int matched_now = count_if(cbegin(genZps), cend(genZps), [](const sidm::Zp& p){return p.matched;});
                    if (matched_now == 2) {
                        zp_r_ = sidm::Zp(q);
                        zp_r_._eventId = eventNum_;
                        darkPhoton_rereco_->Fill();
                        // Store this faking jet
                        j_fake_e_ = sidm::Jet(q.second);
                        jet_properties_->Fill();
                        // Remove the matched pat pair from pat::Electron and pat::Jet collection.
                        sidm::remove_from_collection(&patElectronPtr_, q.first);
                        sidm::remove_from_collection(&patJetPtr_, q.second);

                        PAIR = 2;            //< 2 pairs are reconstructed.
                        EJPAIR = 1;          //< 1 pair is electron-jet.
                        JJPAIR = 0;
                        ++event2pairs_;
                        ++event2pairsEpEj_;
                    }
                }
                if (PAIR == 1)
                    EJPAIR = 0;  //< PAIR remains as same (0/1).
            } else {  //< PAIR==0
                bool first_genZp_matched(false);
                bool second_genZp_matched(false);
                for (const auto& q : ejPairTmp) {
                    sidm::match_patPair_with_zps(q, genZps, dRusb_);
                    if (!first_genZp_matched && genZps[0].matched) {
                        first_genZp_matched = true;
                        zp_r_ = sidm::Zp(q);
                        zp_r_._eventId = eventNum_;
                        darkPhoton_rereco_->Fill();
                        // Store this faking jet
                        j_fake_e_ = sidm::Jet(q.second);
                        jet_properties_->Fill();
                        // Remove the matched pat pair from pat::Electron and pat::Jet collection.
                        sidm::remove_from_collection(&patElectronPtr_, q.first);
                        sidm::remove_from_collection(&patJetPtr_, q.second);
                    }
                    if (!second_genZp_matched && genZps[1].matched) {
                        second_genZp_matched = true;
                        zp_r_ = sidm::Zp(q);
                        zp_r_._eventId = eventNum_;
                        darkPhoton_rereco_->Fill();
                        // Store this faking jet
                        j_fake_e_ = sidm::Jet(q.second);
                        jet_properties_->Fill();
                        // Remove the matched pat pair from pat::Electron and pat::Jet collection.
                        sidm::remove_from_collection(&patElectronPtr_, q.first);
                        sidm::remove_from_collection(&patJetPtr_, q.second);
                    }
                    /// If both matched, skip the rest.
                    if (first_genZp_matched && second_genZp_matched) break;
                }
                int matched_now = static_cast<int>(first_genZp_matched) + static_cast<int>(second_genZp_matched);
                if (matched_now == 2) {
                    PAIR   = 2;          //< 2 Pairs are reconstructed.
                    EJPAIR = 2;          //< Both are electron-jet pairs.
                    JJPAIR = 0;
                    ++event2pairs_;
                    ++event2pairsEjEj_;
                }
                else if (matched_now == 1) { PAIR = 1; EJPAIR = 1; }
                else if (matched_now == 0) { PAIR = 0; EJPAIR = 0; }
            }
        }
    }
    // ---------------EJ SEARCH END---------------

    /**
     * __JJ PAIR SEARCH__
     */
    if (PAIR<2) {
        vector<pair<Ptr<pat::Jet>, Ptr<pat::Jet> > > jjPairTmp{};
        if (patJetPtr_.size() > 1) {  //< it does not make sense if there are less than 2 patJets.
            for (auto it = cbegin(patJetPtr_); it != cend(patJetPtr_)-1; ++it) {
                for(auto jt = cbegin(patJetPtr_)+1; jt != cend(patJetPtr_); ++jt) {
                    jjPairTmp.emplace_back(std::make_pair(*it, *jt));
                }
            }
        }

        // Narrow mass sideband
        jjPairTmp.erase(remove_if(begin(jjPairTmp), end(jjPairTmp), [this](const auto& p){
                    float m = (p.first->p4()+p.second->p4()).M();
                    return m<(1-zpMassSb_)*zpMass_ || m>(1+zpMassSb_)*zpMass_;
                    }), jjPairTmp.end());
        
        if (jjPairTmp.size()==0) {
            // Again no pairs in mass band...
            JJPAIR = 0;
            if (PAIR == 0) ++event0pair_;
            else { // PAIR == 1
                ++event1pair_;
                if (EPPAIR) ++event1pairEp_;
                if (EJPAIR) ++event1pairEj_;
            }
        } else {
            // First, sort by invM
            sort(begin(jjPairTmp), end(jjPairTmp), [this](const auto& lhs, const auto& rhs) {
                    float lhsM = (lhs.first->p4()+lhs.second->p4()).M();
                    float rhsM = (rhs.first->p4()+rhs.second->p4()).M();
                    return std::abs(lhsM-zpMass_) < std::abs(rhsM-zpMass_);
                    });
            sidm::pairvec<pat::Jet> jjPair(jjPairTmp, true);
            jjPairTmp.clear();
            jjPairTmp = jjPair.get_zip();

            if (PAIR == 1) {
                for (const auto& q : jjPairTmp) {
                    sidm::match_patPair_with_zps(jjPairTmp[0], genZps, dRusb_);
                    int matched_now = count_if(cbegin(genZps), cend(genZps), [](const sidm::Zp& p){return p.matched;});
                    if (matched_now == 2) {
                        zp_r_ = sidm::Zp(q);
                        zp_r_._eventId = eventNum_;
                        darkPhoton_rereco_->Fill();
                        // Store these 2 faking jet
                        j_fake_e_ = sidm::Jet(q.first);
                        jet_properties_->Fill();
                        j_fake_e_ = sidm::Jet(q.second);
                        jet_properties_->Fill();
                        // Remove the matched pat pair from pat::Jet collection.
                        sidm::remove_from_collection(&patJetPtr_, q.first);
                        sidm::remove_from_collection(&patJetPtr_, q.second);

                        PAIR = 2;        //< 2 pairs are reconstructed.
                        JJPAIR =1;       //< 1 pair is jet-jet.
                        ++event2pairs_;
                        if (EPPAIR) ++event2pairsEpJj_;
                        if (EJPAIR) ++event2pairsEjJj_;
                    }
                }
                if (PAIR == 1) {
                    JJPAIR = 0;  //< PAIR remains as same (1).
                    ++event1pair_;
                    if (EPPAIR) ++event1pairEp_;
                    if (EJPAIR) ++event1pairEj_;
                }
            } else {  //< PAIR=0
                bool first_genZp_matched(false);
                bool second_genZp_matched(false);
                for (const auto& q : jjPairTmp) {
                    sidm::match_patPair_with_zps(q, genZps, dRusb_);
                    if (!first_genZp_matched && genZps[0].matched) {
                        first_genZp_matched = true;
                        zp_r_ = sidm::Zp(q);
                        zp_r_._eventId = eventNum_;
                        darkPhoton_rereco_->Fill();
                        // Store these 2 faking jet
                        j_fake_e_ = sidm::Jet(q.first);
                        jet_properties_->Fill();
                        j_fake_e_ = sidm::Jet
                        // Remove the matched pat pair from pat::Jet collection.
                        sidm::remove_from_collection(&patJetPtr_, q.first);
                        sidm::remove_from_collection(&patJetPtr_, q.second);
                    }
                    if (!second_genZp_matched && genZps[1].matched) {
                        second_genZp_matched = true;
                        zp_r_ = sidm::Zp(q);
                        zp_r_._eventId = eventNum_;
                        darkPhoton_rereco_->Fill();
                        // Store these 2 faking jet
                        j_fake_e_ = sidm::Jet(q.first);
                        jet_properties_->Fill();
                        j_fake_e_ = sidm::Jet
                        // Remove the matched pat pair from pat::Jet collection.
                        sidm::remove_from_collection(&patJetPtr_, q.first);
                        sidm::remove_from_collection(&patJetPtr_, q.second);
                    }
                    /// If both matched, skip the rest.
                    if (first_genZp_matched && second_genZp_matched) break;
                }
                int matched_now = static_cast<int>(first_genZp_matched) + static_cast<int>(second_genZp_matched);
                if (matched_now == 2) {
                    PAIR = 2;        //< 2 pairs are reconstructed.
                    JJPAIR = 2;      //< Both are jet-jet pairs.
                    ++event2pairs_;
                    ++event2pairsJjJj_;
                } else if (matched_now == 1) {
                    PAIR = 1;        //< 1 pair is reconstructed.
                    JJPAIR = 1;      //< is a jet-jet pair.
                    ++event1pair_;
                    ++event1pairJj_;
                } else if (matched_now == 0) {
                    JJPAIR = 0;  //< PAIR remains as same(0).
                    ++event0pair_;
                }
            }
        }
    }
    // ---------------JJ SEARCH END---------------
    
    assert(PAIR == EPPAIR + EJPAIR + JJPAIR);

    //cout << "[" << eventNum_ << "] "
    //     << "Num.Pair: "<< PAIR << ", EP: " << EPPAIR << ", EJ: " << EJPAIR << ", JJ: " << JJPAIR <<endl;

    ++eventNum_;
}


// ------------ method called once each job just before starting event loop  ------------
void 
sidm::mcValidation::beginJob()
{

    /// Statistics of multiplicities of collections per event
    eventTree_ = fs_->make<TTree>("eventTree", "information per event");
    
    eventTree_->Branch("numberOfElectrons", &electron_N, "numberOfElectrons/I");
    eventTree_->Branch("numberOfPositrons", &positron_N, "numberOfPositrons/I");
    eventTree_->Branch("numberOfElectronsZp", &electron_from_zp_N, "numberOfElectronsZp/I");
    eventTree_->Branch("numberOfPositronsZp", &positron_from_zp_N, "numberOfPositronsZp/I");
    eventTree_->Branch("numberOfZps", &zp_N, "numberOfZps/I");
    eventTree_->Branch("numberOfPs", &ps_N, "numberOfPs/I");
    eventTree_->Branch("numberOfPatElectrons", &patE_N, "numberOfPatElectrons/I");
    eventTree_->Branch("numberOfPatJets", &patJet_N, "numberOfPatJets/I");

    /// darkphoton tree from genParticles
    darkPhoton_reco_ = fs_->make<TTree>("darkPhoton_reco", "gen darkPhotons");

    darkPhoton_reco_->Branch("eventId", &zp_._eventId, "eventId/I");
    darkPhoton_reco_->Branch("pt", &zp_._pt, "pt/F");
    darkPhoton_reco_->Branch("eta", &zp_._eta, "eta/F");
    darkPhoton_reco_->Branch("mass", &zp_._mass, "mass/F");
    darkPhoton_reco_->Branch("invm", &zp_._invM, "invm/F");
    darkPhoton_reco_->Branch("pt_e", &zp_.e._pt, "pt_e/F");
    darkPhoton_reco_->Branch("pt_p", &zp_.p._pt, "pt_p/F");
    darkPhoton_reco_->Branch("dR_ep", &zp_._dR, "dR_ep/F");
    darkPhoton_reco_->Branch("dEta_ep", &zp_._dEta, "dEta_ep/F");
    darkPhoton_reco_->Branch("dPhi_ep", &zp_._dPhi, "dPhi_ep/F");
    darkPhoton_reco_->Branch("dv_x", &zp_._dv_x, "dv_x/F");
    darkPhoton_reco_->Branch("dv_y", &zp_._dv_y, "dv_y/F");
    darkPhoton_reco_->Branch("dv_z", &zp_._dv_z, "dv_z/F");

    /// pseudo-scalor tree from genParticles
    pscalar_reco_ = fs_->make<TTree>("pscalar_reco", "gen pseudo-scalars");

    pscalar_reco_->Branch("eventId", &ps_._eventId, "eventId/I");
    pscalar_reco_->Branch("mass", &ps_._mass, "mass/F");
    pscalar_reco_->Branch("invm", &ps_._invM, "invm/F");
    pscalar_reco_->Branch("dEta", &ps_._dEta, "dEta/F");
    pscalar_reco_->Branch("dPhi", &ps_._dPhi, "dPhi/F");
    pscalar_reco_->Branch("dR", &ps_._dR, "dR/F");

    /// darkphoton tree from reco particles
    darkPhoton_rereco_ = fs_->make<TTree>("darkPhoton_rereco", "reco darkPhotons");

    darkPhoton_rereco_->Branch("eventId", &zp_r_._eventId, "eventId/I");
    darkPhoton_rereco_->Branch("pt", &zp_r_._pt, "pt/F");
    darkPhoton_rereco_->Branch("eta", &zp_r_._eta, "eta/F");
    darkPhoton_rereco_->Branch("invm", &zp_r_._invM, "invm/F");
    darkPhoton_rereco_->Branch("pt_e", &zp_r_.e._pt, "pt_e/F");
    darkPhoton_rereco_->Branch("pt_p", &zp_r_.p._pt, "pt_p/F");
    darkPhoton_rereco_->Branch("dR_ep", &zp_r_._dR, "dR_ep/F");
    darkPhoton_rereco_->Branch("dEta_ep", &zp_r_._dEta, "dEta_ep/F");
    darkPhoton_rereco_->Branch("dPhi_ep", &zp_r_._dPhi, "dPhi_ep/F");

    /// jets that could be mis-reconstructed electrons
    jet_properties_ = fs_->make<TTree>("jet_properties", "mis-recoed jets properties");

    jet_properties_->Branch("eventId", &j_fake_e_._eventId, "eventId/I");
    jet_properties_->Branch("pt", &j_fake_e_._pt, "pt/F");
    jet_properties_->Branch("eta", &j_fake_e_._eta, "eta/F");
    jet_properties_->Branch("charged_H_over_E", &j_fake_e_._charged_H_over_E, "charged_H_over_E/F");
    jet_properties_->Branch("H_over_E", &j_fake_e_._H_over_E, "H_over_E/F");
    jet_properties_->Branch("electronEnergyFraction", &j_fake_e_._electronEnergyFraction, "electronEnergyFraction/F");
    jet_properties_->Branch("chargedMultiplicity", &j_fake_e_._chargedMultiplicity, "chargedMultiplicity/I");
    jet_properties_->Branch("electronMultiplicity", &j_fake_e_._electronMultiplicity, "electronMultiplicity/I");
    jet_properties_->Branch("num_elec_fullin", &j_fake_e_._num_elec_fullin, "num_elec_fullin/I");
    jet_properties_->Branch("num_elec_curledin", &j_fake_e_._num_elec_curledin, "num_elec_curledin/I");
    jet_properties_->Branch("num_elec_curledout", &j_fake_e_._num_elec_curledout, "num_elec_curledout/I");    
}

// ------------ method called once each job just after ending the event loop  ------------
void 
sidm::mcValidation::endJob() 
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
sidm::mcValidation::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(sidm::mcValidation);
