// -*- C++ -*-
//
// Package:    sidm-bs/mcValidation
// Class:      mcValidation
#include <algorithm>
#include <cmath>

#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/Math/interface/LorentzVector.h"
#include "DataFormats/Math/interface/deltaR.h"

#include "sidm-bs/jetComponent/interface/physicsObject.h"
#include "sidm-bs/jetComponent/interface/mcValidation.h"

sidm::mcValidation::mcValidation(const edm::ParameterSet& iConfig):
    genParticleTk_(consumes<edm::View<reco::GenParticle> >(iConfig.getUntrackedParameter<edm::InputTag>("GenParticleTag_", edm::InputTag("prunedGenParticles")))),
    ssVerticeTk_(consumes<edm::View<reco::VertexCompositePtrCandidate> >(iConfig.getUntrackedParameter<edm::InputTag>("SsVerticeTag_", edm::InputTag("slimmedSecondaryVertices")))),
    patElectronTk_(consumes<edm::View<pat::Electron> >(iConfig.getUntrackedParameter<edm::InputTag>("PatElectronTag_", edm::InputTag("slimmedElectrons")))),
    patJetTk_(consumes<edm::View<pat::Jet> >(iConfig.getUntrackedParameter<edm::InputTag>("PatJetTag_", edm::InputTag("slimmedJets")))),
    zpMassSb_(iConfig.getUntrackedParameter<double>("ZpMassSideBand"))
{
    usesResource("TFileService");
    eventNum_ = 0;
}


sidm::mcValidation::~mcValidation()
{

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

    vector<pair<const reco::Candidate*, const reco::Candidate*> > genZps{}; // {<e,p>, <e,p>}
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

                genZps.emplace_back(make_pair(p.daughter(0), p.daughter(1)));

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

                genZps.emplace_back(make_pair(p.daughter(1), p.daughter(0)));
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
            ps_._dEta    =  std::abs( p.daughter(0)->eta() - p.daughter(1)->eta() );
            ps_._dPhi    =  std::abs( p.daughter(0)->phi() - p.daughter(1)->phi() );
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

    vector<pat::Electron> patEp_es;
    vector<pat::Electron> patEp_ps;
    copy_if(patElectronHdl_->begin(), patElectronHdl_->end(), back_inserter(patEp_es), [](const pat::Electron p){return p.charge()<0;});
    copy_if(patElectronHdl_->begin(), patElectronHdl_->end(), back_inserter(patEp_ps), [](const pat::Electron p){return p.charge()>0;});
    //assert(patEp_es.size()>0 && patEp_ps.size()>0);
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
        }
    } else {
        cout << "[" << eventNum_ << "] e:"
             << patEp_es.size() << " p: " << patEp_ps.size() << endl;
    }

    ++eventNum_;

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
