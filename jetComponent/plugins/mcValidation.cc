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
}


sidm::mcValidation::~mcValidation()
{

    eventNum_ = 0;
    //sidm::test ss(88);
    //ss.get();

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

    electron_N = std::count_if(cbegin(*genParticleHdl_),
                               cend(*genParticleHdl_),
                               [](const auto& p){ return p.pdgId() == 11 && 
                                                         p.status() == 1 &&
                                                         p.pt() > 2. ; });
    electron_from_zp_N = std::count_if(cbegin(*genParticleHdl_),
                                       cend(*genParticleHdl_),
                                       [](const auto& p){ return p.pdgId() == 11 &&
                                                                 p.mother()->pdgId() == 32; });

    positron_N = std::count_if(cbegin(*genParticleHdl_),
                               cend(*genParticleHdl_),
                               [](const auto& p){ return p.pdgId() == -11 &&
                                                         p.status() == 1  &&
                                                         p.pt() > 2. ; });
    positron_from_zp_N = std::count_if(cbegin(*genParticleHdl_),
                                       cend(*genParticleHdl_),
                                       [](const auto& p){ return p.pdgId() == -11 &&
                                                                 p.mother()->pdgId() == 32; });

    zp_N = std::count_if(cbegin(*genParticleHdl_),
                         cend(*genParticleHdl_),
                         [](const auto& p){ return p.pdgId() == 32; });
    ps_N = std::count_if(cbegin(*genParticleHdl_),
                         cend(*genParticleHdl_),
                         [](const auto& p){ return p.pdgId() == 35; });

    patE_N = patElectronHdl_->size();
    patJet_N = patJetHdl_->size();

    eventTree_->Fill();

    //------------------------------------------
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

            if ( p.daughter(0)->pdgId() *
                    p.daughter(1)->pdgId() != -11*11 ) {
                cout << "[" << eventNum_ << "] "
                     << "Found a darkphoton whose daughters are not electron pair" << endl;
                cout << "--> Instead (" << p.daughter(0)->pdgId() << ", "
                                        << p.daughter(1)->pdgId() << ")" <<endl;

                continue;
            }

            if ( p.daughter(0)->pdgId() == 11 ) {
                // 0: electron, 1: positron
                zp_.e._eta     = p.daughter(0)->eta();
                zp_.e._phi     = p.daughter(0)->phi();
                zp_.e._pt      = p.daughter(0)->pt();
                zp_.e._energy  = p.daughter(0)->energy();

                zp_.p._eta     = p.daughter(1)->eta();
                zp_.p._phi     = p.daughter(1)->phi();
                zp_.p._pt      = p.daughter(1)->pt();
                zp_.p._energy  = p.daughter(1)->energy();

                genZps.push_back(zp_);

            } else {
                // 0: positron, 1: electron
                zp_.p._eta     = p.daughter(0)->eta();
                zp_.p._phi     = p.daughter(0)->phi();
                zp_.p._pt      = p.daughter(0)->pt();
                zp_.p._energy  = p.daughter(0)->energy();

                zp_.e._eta     = p.daughter(1)->eta();
                zp_.e._phi     = p.daughter(1)->phi();
                zp_.e._pt      = p.daughter(1)->pt();
                zp_.e._energy  = p.daughter(1)->energy();

                genZps.push_back(zp_);
            }

            zp_._eventId = eventNum_;
            zp_._pt      = p.pt();
            zp_._eta     = p.eta();
            zp_._mass    = p.mass();
            zp_._invM    = ( p.daughter(0)->p4() + p.daughter(1)->p4() ).M();
            zp_._dR      = deltaR( *(p.daughter(0)), *(p.daughter(1)) );
            zp_._dEta    = std::abs( p.daughter(0)->eta() - p.daughter(1)->eta() );
            zp_._dPhi    = std::abs( p.daughter(0)->phi() - p.daughter(1)->phi() );
            zp_._dv_x    = p.daughter(0)->vx();
            zp_._dv_y    = p.daughter(0)->vy();
            zp_._dv_z    = p.daughter(0)->vz();

            darkPhoton_reco_->Fill(); 

        }

        /// Ps
        if ( p.pdgId() == 35 ) {
        
            if ( p.numberOfDaughters() != 2 ) { 
                continue;
            }

            if ( p.daughter(0)->pdgId() *
                 p.daughter(1)->pdgId() != 32*32 ) {
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

    int PAIR(-1), EPPAIR(-1), EJPAIR(-1), JJPAIR(-1);

    vector<Ptr<pat::Electron> > patElectronPtr_ = patElectronHdl_->ptrs();
    vector<Ptr<pat::Electron> > patEp_es{}, patEp_ps{};

    copy_if(patElectronPtr_.begin(), patElectronPtr_.end(), back_inserter(patEp_es), [](const Ptr<pat::Electron>& p){return p->charge()<0;});
    copy_if(patElectronPtr_.begin(), patElectronPtr_.end(), back_inserter(patEp_ps), [](const Ptr<pat::Electron>& p){return p->charge()>0;});
    //assert(patEp_es.size()>0 && patEp_ps.size()>0);

    /* EP PAIR SEARCH */
    if ( patEp_es.size()>0 && patEp_ps.size()>0 ) {

        // construct pair vector
        pairvec<pat::Electron> epPair(patEp_es, patEp_ps);
        vector<pair<Ptr<pat::Electron>, Ptr<pat::Electron> > > epPairTmp = epPair.get();

        // narrow mass sideband
        remove_if(begin(epPairTmp), end(epPairTmp), [this](const auto& p){ 
                float m = (p.first->p4()+p.second->p4()).M();
                return m>=(1-zpMassSb_)*zpMass_ && m<=(1+zpMassSb_)*zpMass_;
                });

        if ( epPairTmp.size() == 0 ) {
            PAIR   = 0;
            EPPAIR = 0;
        } else if (epPairTmp.size()==1) {
            //ONLY 1 pair in mass band 
            //--> looking at dR matched or not with one of the two dark photons.
            vector<float> avedR{};
            for (auto& p : genZps) {
                std::pair<float, float> drTmp(p.dRVal(epPairTmp[0]));
                if (drTmp.first <= dRusb_ && drTmp.second <= dRusb_)
                    p.matched = true;
                avedR.push_back((drTmp.first+drTmp.second)/2);
            }
            
            int matchedNum = count_if(cbegin(genZps), cend(genZps), [](const sidm::Zp& p){return p.matched;});

            if (matchedNum == 2) {
               //whose average dR is smaller got matched, the other one remains as unmatched.
               if (avedR[0] <= avedR[1])
                   genZps[1].matched = false;
               else
                   genZps[0].matched = false;
               matchedNum = 1;
            }
            if (matchedNum == 1) {
                PAIR   = 1; EPPAIR = 1;

                zp_r_ = sidm::Zp(epPairTmp[0]);
                zp_r_._eventId = eventNum_;
                darkPhoton_rereco_->Fill();
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
            pairvec<pat::Electron> epPair(epPairTmp, true);
            epPairTmp.clear();
            epPairTmp = epPair.get_zip();

            if (epPairTmp.size() == 1) {
                //ONLY 1 unique pair, esscencially the same way as we deal with 1 pair
                vector<float> avedR{};
                for (auto& p : genZps) {
                    std::pair<float, float> drTmp(p.dRVal(epPairTmp[0]));
                    if (drTmp.first <= dRusb_ && drTmp.second <= dRusb_)
                        p.matched = true;
                    avedR.push_back((drTmp.first+drTmp.second)/2);
                }

                int matchedNum = count_if(cbegin(genZps), cend(genZps), [](const sidm::Zp& p){return p.matched;});

                if (matchedNum == 2) {
                    //whose average dR is larger got matched, the other one remains as unmatched.
                    if (avedR[0] <= avedR[1])
                        genZps[1].matched = false;
                    else
                        genZps[0].matched = false;
                    matchedNum = 1;
                }

                if (matchedNum == 1) {
                    PAIR = 1; EPPAIR = 1;
                    zp_r_ = sidm::Zp(epPairTmp[0]);
                    zp_r_._eventId = eventNum_;
                    darkPhoton_rereco_->Fill();
                } else { PAIR = 0; EPPAIR = 0; }

            } else {
                //More than 1 unique pair in mass band-
                for (const auto& q : epPairTmp) {
                    if (genZps[0].matched && genZps[1].matched) break;
                    std::pair<float, float> drTmp0(genZps[0].dRVal(q));
                    std::pair<float, float> drTmp1(genZps[1].dRVal(q));
                    if (!genZps[0].matched && drTmp0.first <= dRusb_ && drTmp0.second <= dRusb_) {
                        genZps[0].matched = true;
                        zp_r_ = sidm::Zp(q);
                        zp_r_._eventId = eventNum_;
                        darkPhoton_rereco_->Fill();
                        continue;
                    }
                    if (!genZps[1].matched && drTmp1.first <= dRusb_ && drTmp1.second <= dRusb_) {
                        genZps[1].matched = true;
                        zp_r_ = sidm::Zp(q);
                        zp_r_._eventId = eventNum_;
                        darkPhoton_rereco_->Fill();
                    }
                }

                int matchedNum = count_if(cbegin(genZps), cend(genZps), [](const sidm::Zp& p){return p.matched;});
                if (matchedNum == 2) {
                    PAIR = 2;
                    EPPAIR = 2;
                    ++event2pairs_;
                    ++event2pairsEpEp_;
                }
                else if (matchedNum == 1) { PAIR = 1; EPPAIR = 1; }
                else if (matchedNum == 0) { PAIR = 0; EPPAIR = 0; }
                
            }
        }
            
    } else if ( patEp_es.size() >0 || patEp_ps.size()>0 ) {
    
    } else {
    
    }
    /*--
    if (patEp_es.size()>=2 && patEp_ps.size()>=2) {
        for (auto& p : genZps) {
            // sort by dR separately for electrons & positrons, from large to small 
            sort(begin(patEp_es), end(patEp_es), [p](const pat::Electron& lhs, const pat::Electron& rhs){return deltaR(lhs, *p.first) > deltaR(rhs, *p.first);});
            sort(begin(patEp_ps), end(patEp_ps), [p](const pat::Electron& lhs, const pat::Electron& rhs){return deltaR(lhs, *p.second) > deltaR(rhs, *p.second);});
            
            zp_r_._eventId  = eventNum_;

            zp_r_.e._eta    = patEp_es[-1].eta();
            zp_r_.e._phi    = patEp_es[-1].phi();
            zp_r_.e._pt     = patEp_es[-1].pt();
            zp_r_.e._energy = patEp_es[-1].energy();
            zp_r_.p._eta    = patEp_ps[-1].eta();
            zp_r_.p._phi    = patEp_ps[-1].phi();
            zp_r_.p._pt     = patEp_ps[-1].pt();
            zp_r_.p._energy = patEp_ps[-1].energy();

            zp_r_._eta  = zp_r_.e._eta + zp_r_.p._eta;
            zp_r_._invM = ( patEp_es[-1].p4(), patEp_ps[-1].p4() ).M();
            zp_r_._dEta = abs( zp_r_.e._eta - zp_r_.p._eta );
            zp_r_._dPhi = abs( zp_r_.p._phi - zp_r_.p._phi );
            zp_r_._dR   = deltaR( patEp_es[-1], patEp_ps[-1] );

            darkPhoton_rereco_->Fill();

            patEp_es.pop_back();
            patEp_ps.pop_back();

            cout << "[" << eventNum_ << "] invm:"
                 << zp_r_._invM << endl;
        }
    } else {
        //cout << "[" << eventNum_ << "] e:"
        //     << patEp_es.size() << " p: " << patEp_ps.size() << endl;
    }
    --*/

    ++eventNum_;
    cout << PAIR << EPPAIR << EJPAIR << JJPAIR <<endl;

}


// ------------ method called once each job just before starting event loop  ------------
void 
sidm::mcValidation::beginJob()
{

    eventTree_ = fs_->make<TTree>("eventTree", "information per event");
    eventTree_->Branch("numberOfElectrons", &electron_N, "numberOfElectrons/I");
    eventTree_->Branch("numberOfPositrons", &positron_N, "numberOfPositrons/I");
    eventTree_->Branch("numberOfElectronsZp", &electron_from_zp_N, "numberOfElectronsZp/I");
    eventTree_->Branch("numberOfPositronsZp", &positron_from_zp_N, "numberOfPositronsZp/I");
    eventTree_->Branch("numberOfZps", &zp_N, "numberOfZps/I");
    eventTree_->Branch("numberOfPs", &ps_N, "numberOfPs/I");
    eventTree_->Branch("numberOfPatElectrons", &patE_N, "numberOfPatElectrons/I");
    eventTree_->Branch("numberOfPatJets", &patJet_N, "numberOfPatJets/I");

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

    pscalar_reco_ = fs_->make<TTree>("pscalar_reco", "gen pseudo-scalars");

    pscalar_reco_->Branch("eventId", &ps_._eventId, "eventId/I");
    pscalar_reco_->Branch("mass", &ps_._mass, "mass/F");
    pscalar_reco_->Branch("invm", &ps_._invM, "invm/F");
    pscalar_reco_->Branch("dEta", &ps_._dEta, "dEta/F");
    pscalar_reco_->Branch("dPhi", &ps_._dPhi, "dPhi/F");
    pscalar_reco_->Branch("dR", &ps_._dR, "dR/F");

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
