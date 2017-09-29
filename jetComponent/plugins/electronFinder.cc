// -*- C++ -*-
//
// Package:    sidm-bs/electronFinder
// Class:      electronFinder
#include <algorithm>
#include <cmath>

#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/Common/interface/Ptr.h"
#include "DataFormats/Math/interface/LorentzVector.h"
#include "DataFormats/Math/interface/deltaR.h"

#include "sidm-bs/jetComponent/interface/physicsObject.h"
#include "sidm-bs/jetComponent/interface/utilities.h"

#include "sidm-bs/jetComponent/interface/electronFinder.h"

sidm::electronFinder::electronFinder(const edm::ParameterSet& iConfig):
    genParticleTk_(consumes<edm::View<reco::GenParticle> >(iConfig.getUntrackedParameter<edm::InputTag>("GenParticleTag_", edm::InputTag("prunedGenParticles")))),
    pfTk_(consumes<edm::View<pat::PackedCandidate> >(iConfig.getUntrackedParameter<edm::InputTag>("PatPfTag_", edm::InputTag("packedPFCandidates")))),
    patElectronTk_(consumes<edm::View<pat::Electron> >(iConfig.getUntrackedParameter<edm::InputTag>("PatElectronTag_", edm::InputTag("slimmedElectrons")))),
    patJetTk_(consumes<edm::View<pat::Jet> >(iConfig.getUntrackedParameter<edm::InputTag>("PatJetTag_", edm::InputTag("slimmedJets"))))
{
    usesResource("TFileService");
    eventNum_ = 0;
}


sidm::electronFinder::~electronFinder()
{

    std::cout << "\n\nNumber of events: " << eventNum_ << "\n\n";
    eventNum_ = 0;

}

// ------------ method called for each event  ------------
void
sidm::electronFinder::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
    using namespace edm;
    using namespace std;

    Handle<View<reco::GenParticle> > genParticleHdl_;
    iEvent.getByToken(genParticleTk_, genParticleHdl_);

    Handle<View<pat::Electron> > patElectronHdl_;
    iEvent.getByToken(patElectronTk_, patElectronHdl_);

    Handle<View<pat::Jet> > patJetHdl_;
    iEvent.getByToken(patJetTk_, patJetHdl_);

    Handle<View<pat::PackedCandidate> > pfHdl_;
    iEvent.getByToken(pfTk_, pfHdl_);

    pfElectron_N = count_if(cbegin(*pfHdl_), cend(*pfHdl_),
            [](const auto& p){ return abs(p.pdgId()) == 11; });
    //cout<<"Event"<<eventNum_<<": "<<pfElectron_N<<endl;
    eventTree_->Fill();

    ++eventNum_;
}


// ------------ method called once each job just before starting event loop  ------------
void 
sidm::electronFinder::beginJob()
{
    /// Statistics of multiplicities of collections per event
    eventTree_ = fs_->make<TTree>("eventTree", "information per event");
    
    eventTree_->Branch("numOfEleInPfCands", &pfElectron_N, "numOfEleInPfCands/I");
    // eventTree_->Branch("numberOfElectrons", &electron_N, "numberOfElectrons/I");

}

// ------------ method called once each job just after ending the event loop  ------------
void 
sidm::electronFinder::endJob() 
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
sidm::electronFinder::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(sidm::electronFinder);
