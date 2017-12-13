// -*- C++ -*-
//
// Package:    sidm-bs/MuonAnalysis
// Class:      MuonAnalysis
#include <algorithm>
#include <cmath>

#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/Common/interface/Ptr.h"
#include "DataFormats/Math/interface/LorentzVector.h"
#include "DataFormats/Math/interface/deltaR.h"

#include "sidm-bs/jetComponent/interface/physicsObject.h"
#include "sidm-bs/jetComponent/interface/utilities.h"

#include "sidm-bs/jetComponent/interface/MuonAnalysis.h"

sidm::MuonAnalysis::MuonAnalysis(const edm::ParameterSet& iConfig):
    genParticleTk_(consumes<edm::View<reco::GenParticle> >(iConfig.getUntrackedParameter<edm::InputTag>("GenParticleTag_", edm::InputTag("prunedGenParticles")))),
    pfTk_(consumes<edm::View<pat::PackedCandidate> >(iConfig.getUntrackedParameter<edm::InputTag>("PatPfTag_", edm::InputTag("packedPFCandidates")))),
    patElectronTk_(consumes<edm::View<pat::Electron> >(iConfig.getUntrackedParameter<edm::InputTag>("PatElectronTag_", edm::InputTag("slimmedElectrons")))),
    patMuonTk_(consumes<edm::View<pat::Muon> >(iConfig.getUntrackedParameter<edm::InputTag>("PatMuonTag_", edm::InputTag("slimmedMuons")))),
    patJetTk_(consumes<edm::View<pat::Jet> >(iConfig.getUntrackedParameter<edm::InputTag>("PatJetTag_", edm::InputTag("slimmedJets"))))
{
    usesResource("TFileService");
    eventNum_ = 0;
}


sidm::MuonAnalysis::~MuonAnalysis()
{

    //std::cout << "\n\nNumber of events: " << eventNum_ << "\n\n";
    edm::LogDebug("JobDone")<<"Total events: "<<eventNum_;
    eventNum_ = 0;

}

// ------------ method called for each event  ------------
void
sidm::MuonAnalysis::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
    using namespace edm;
    using namespace std;

    Handle<View<reco::GenParticle> > genParticleHdl_;
    iEvent.getByToken(genParticleTk_, genParticleHdl_);

    Handle<View<pat::Electron> > patElectronHdl_;
    iEvent.getByToken(patElectronTk_, patElectronHdl_);

    Handle<View<pat::Muon> > patMuonHdl_;
    iEvent.getByToken(patMuonTk_, patMuonHdl_);

    Handle<View<pat::Jet> > patJetHdl_;
    iEvent.getByToken(patJetTk_, patJetHdl_);

    Handle<View<pat::PackedCandidate> > pfHdl_;
    iEvent.getByToken(pfTk_, pfHdl_);

    //vector<Ptr<pat::Muon> > patMuonPtr_ = patMuonHdl_->ptrs();

    muonMultiplicity_->Fill(patMuonHdl_->size());
    for (const auto& _mu : *patMuonHdl_) {
        muonEnergy_->Fill(_mu.energy());
        muonPt_->Fill(_mu.pt());
        muonEt_->Fill(_mu.et());
    }

    ++eventNum_;
}


// ------------ method called once each job just before starting event loop  ------------
void 
sidm::MuonAnalysis::beginJob()
{
    /// Statistics of multiplicities of collections per event
    eventTree_ = fs_->make<TTree>("eventTree", "information per event");
    
    TFileDirectory muonDir = fs_->mkdir("MuonAlone");
    muonMultiplicity_ = muonDir.make<TH1F>("MuonMultiplicity", "Muon Multiplicity", 10,0,10);
    muonEnergy_       = muonDir.make<TH1F>("MuonEnergy", "Muon Energy", 100,0,100);
    muonPt_           = muonDir.make<TH1F>("MuonPt", "Muon Pt", 100,0,100);
    muonEt_           = muonDir.make<TH1F>("MuonEt", "Muon Et", 100,0,100);
    // eventTree_->Branch("numberOfElectrons", &electron_N, "numberOfElectrons/I");

}

// ------------ method called once each job just after ending the event loop  ------------
void 
sidm::MuonAnalysis::endJob() 
{
    muonMultiplicity_->GetXaxis()->SetTitle("Number of Muons");
    muonMultiplicity_->GetYaxis()->SetTitle("Event Number");
    
    muonEnergy_->GetXaxis()->SetTitle("Energy [GeV]");
    muonEnergy_->GetYaxis()->SetTitle("Event Number");

    muonPt_->GetXaxis()->SetTitle("pT [GeV]");
    muonPt_->GetYaxis()->SetTitle("Event Number");
    
    muonEt_->GetXaxis()->SetTitle("ET [GeV]");
    muonEt_->GetYaxis()->SetTitle("Event Number");

}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
sidm::MuonAnalysis::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(sidm::MuonAnalysis);
