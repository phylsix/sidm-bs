// -*- C++ -*-
//
// Package:    sidm-bs/mcValidation
// Class:      mcValidation
#include <algorithm>

#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/Math/interface/LorentzVector.h"
#include "DataFormats/Math/interface/deltaR.h"

#include "sidm-bs/jetComponent/interface/physicsObject.h"
#include "sidm-bs/jetComponent/interface/mcValidation.h"

sidm::mcValidation::mcValidation(const edm::ParameterSet& iConfig):
    genParticleTk_(consumes<edm::View<reco::GenParticle> >(iConfig.getUntrackedParameter<edm::InputTag>("GenParticleTag_", edm::InputTag("prunedGenParticles"))))
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

    electron_N = std::count_if(cbegin(*genParticleHdl_),
                               cend(*genParticleHdl_),
                               [](const auto& p){ return p.pdgId() == 11; });
    electron_from_zp_N = std::count_if(cbegin(*genParticleHdl_),
                                       cend(*genParticleHdl_),
                                       [](const auto& p){ return p.pdgId() == 11 &&
                                                                 p.mother()->pdgId() == 32; });

    positron_N = std::count_if(cbegin(*genParticleHdl_),
                               cend(*genParticleHdl_),
                               [](const auto& p){ return p.pdgId() == -11; });
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

    /*--
    cout << "[" << eventNum_ << "] Status (";
    for_each(cbegin(*genParticleHdl_),
             cend(*genParticleHdl_),
             [](const auto& p){
              if (abs(p.pdgId()) == 11 && p.mother()->pdgId() != 32)
                cout << p.status() << ", ";});
    cout << ")" << endl;
    --*/
    eventTree_->Fill();

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

            }

            zp_._eventId = eventNum_;
            zp_._pt = p.pt();
            zp_._mass = p.mass();
            zp_._invM = ( p.daughter(0)->p4() +
                          p.daughter(1)->p4() ).M();
            zp_._dR_ep = deltaR( *(p.daughter(0)),
                                 *(p.daughter(1)) );

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
            ps_._mass = p.mass();
            ps_._invM = ( p.daughter(0)->p4() +
                          p.daughter(1)->p4() ).M();
            pscalar_reco_->Fill();

        }
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

    darkPhoton_reco_ = fs_->make<TTree>("darkPhoton_reco", "gen darkPhotons");

    darkPhoton_reco_->Branch("eventId", &zp_._eventId, "eventId/I");
    darkPhoton_reco_->Branch("pt", &zp_._pt, "pt/F");
    darkPhoton_reco_->Branch("mass", &zp_._mass, "mass/F");
    darkPhoton_reco_->Branch("invm", &zp_._invM, "invm/F");
    darkPhoton_reco_->Branch("pt_e", &zp_.e._pt, "pt_e/F");
    darkPhoton_reco_->Branch("pt_p", &zp_.p._pt, "pt_p/F");
    darkPhoton_reco_->Branch("dR_ep", &zp_._dR_ep, "dR_ep/F");

    pscalar_reco_ = fs_->make<TTree>("pscalar_reco", "gen pseudo-scalars");
    pscalar_reco_->Branch("eventId", &ps_._eventId, "eventId/I");
    pscalar_reco_->Branch("mass", &ps_._mass, "mass/F");
    pscalar_reco_->Branch("invm", &ps_._invM, "invm/F");

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
