// -*- C++ -*-
//
// Package:    sidm-bs/MuonAnalysis
// Class:      MuonAnalysis
#include <algorithm>
#include <cmath>
#include <map>

#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/Common/interface/Ptr.h"
#include "DataFormats/Math/interface/LorentzVector.h"
#include "DataFormats/Math/interface/deltaR.h"

#include "sidm-bs/jetComponent/interface/physicsObject.h"
#include "sidm-bs/jetComponent/interface/utilities.h"

#include "sidm-bs/jetComponent/interface/MuonAnalysis.h"

sidm::MuonAnalysis::MuonAnalysis(const edm::ParameterSet& iConfig):
    genParticleTk_(consumes<edm::View<reco::GenParticle> >(iConfig.getUntrackedParameter<edm::InputTag>("GenParticleTag_", edm::InputTag("prunedGenParticles")))),
    pkdGenTk_(consumes<edm::View<pat::PackedGenParticle> >(iConfig.getUntrackedParameter<edm::InputTag>("PatPkdGenTag_", edm::InputTag("packedGenParticles")))),
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

    std::cout << "\n\nNumber of events: " << eventNum_ << "\n\n";
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

    Handle<View<pat::PackedGenParticle> > pkdGenHdl_;
    iEvent.getByToken(pkdGenTk_, pkdGenHdl_);

    Handle<View<pat::Electron> > patElectronHdl_;
    iEvent.getByToken(patElectronTk_, patElectronHdl_);

    Handle<View<pat::Muon> > patMuonHdl_;
    iEvent.getByToken(patMuonTk_, patMuonHdl_);

    Handle<View<pat::Jet> > patJetHdl_;
    iEvent.getByToken(patJetTk_, patJetHdl_);

    Handle<View<pat::PackedCandidate> > pfHdl_;
    iEvent.getByToken(pfTk_, pfHdl_);

    //vector<Ptr<pat::Muon> > patMuonPtr_ = patMuonHdl_->ptrs();
    vector<Ptr<reco::GenParticle> >      genPtr_    = genParticleHdl_->ptrs();
    vector<Ptr<pat::PackedGenParticle> > pkdGenPtr_ = pkdGenHdl_->ptrs();


    muonMultiplicity_->Fill(patMuonHdl_->size());
    for (const auto& _mu : *patMuonHdl_) {
        muonEnergy_->Fill(_mu.energy());
        muonPt_->Fill(_mu.pt());
        muonEt_->Fill(_mu.et());
        muonEta_->Fill(_mu.eta());
        muonPhi_->Fill(_mu.phi());
        muonVertex_->Fill(sqrt(_mu.vx()*_mu.vx() + _mu.vy()*_mu.vy()), abs(_mu.vz()));
    }

    /**
     * _genDarkphotos_
     * KEY:   The key of darkPhoton in reco::GenParticle collection
     * VALUE: Pair of keys of corresponding muon daughters in pat::PackedGenParticles collection,
     *        first key is muon+, second is muon-.
     * 
     * Keys can be used to uniquely construct edm::Ptr together with Handle so we
     * can keep track.
     */
    map<int, pair<int, int> > genDarkphotons;
    genDarkphotons.clear();

    //cout<<logLabel<<"["<<count_if(cbegin(pkdGenPtr_), cend(pkdGenPtr_), [](const auto& p){return abs(p->pdgId())==13 && p->status()==1;})<<"]\n";
    for (const auto& _dp : genPtr_) {
        if (_dp->pdgId() != 32) continue;
        vector<Ptr<pat::PackedGenParticle> > legitMuonDaughters{};
        
        for (const auto& _gMu : pkdGenPtr_) {
            if (_gMu->status()!=1 || abs(_gMu->pdgId())!=13) continue;
            const reco::Candidate* motherInPrunedCollection(_gMu->mother(0));
            if (motherInPrunedCollection!=nullptr && sidm::is_ancestor(_dp.get(),motherInPrunedCollection)) {
                legitMuonDaughters.push_back(_gMu);
            }
        }

        if (legitMuonDaughters.size()!=2) {
            cout<<"Event"<<setw(5)<<left<<eventNum_<<__LINE__<<"Dark photon has "<<legitMuonDaughters.size()<<" muon daughters instead of 2.\n";
            continue;
        }

        if (legitMuonDaughters[0]->charge()*legitMuonDaughters[1]->charge() != -1) {
            cout<<"Event"<<setw(5)<<left<<eventNum_<<__LINE__<<"Dark photon has same sign muon daughters!\n";
            continue;
        }

        if (legitMuonDaughters[0]->charge() > 0) {
           genDarkphotons[_dp.key()] = std::make_pair(legitMuonDaughters[0].key(), legitMuonDaughters[1].key());
        } else {
           genDarkphotons[_dp.key()] = std::make_pair(legitMuonDaughters[1].key(), legitMuonDaughters[0].key());
        }

        for (const auto& _gMu : legitMuonDaughters) {
            genMuonEnergy_->Fill(_gMu->energy());
            genMuonPt_    ->Fill(_gMu->pt());
            genMuonEt_    ->Fill(_gMu->et());
            genMuonEta_   ->Fill(_gMu->eta());
            genMuonPhi_   ->Fill(_gMu->phi());
            genMuonVertex_->Fill(sqrt(_dp->daughter(0)->vx()*_dp->daughter(0)->vx() +
                                      _dp->daughter(0)->vy()*_dp->daughter(0)->vy()),
                                  abs(_dp->daughter(0)->vz()));
            genMuonVertex_->Fill(sqrt(_dp->daughter(0)->vx()*_dp->daughter(0)->vx() +
                                      _dp->daughter(0)->vy()*_dp->daughter(0)->vy()),
                                  abs(_dp->daughter(0)->vz()));

        }

    }
    genMuonMultiplicity_->Fill(2*genDarkphotons.size());

    if (genDarkphotons.size() != 2) {
        cout<<"Event"<<setw(5)<<left<<eventNum_<<__LINE__<<genDarkphotons.size()<<" gen dark photons found instead of 2.\n";
        return;
    }
    ++eventNum_;
}


// ------------ method called once each job just before starting event loop  ------------
void 
sidm::MuonAnalysis::beginJob()
{
    /// Statistics of multiplicities of collections per event
    eventTree_ = fs_->make<TTree>("eventTree", "information per event");
    
    TFileDirectory muonDir = fs_->mkdir("PAT_Muon");
    muonMultiplicity_ = muonDir.make<TH1F>("MuonMultiplicity", "Muon Multiplicity", 10,0,10);
    muonEnergy_       = muonDir.make<TH1F>("MuonEnergy", "Muon Energy", 100,0,100);
    muonPt_           = muonDir.make<TH1F>("MuonPt", "Muon Pt", 100,0,100);
    muonEt_           = muonDir.make<TH1F>("MuonEt", "Muon Et", 100,0,100);
    muonEta_          = muonDir.make<TH1F>("MuonEta", "Muon Eta",60,-3,3);
    muonPhi_          = muonDir.make<TH1F>("MuonPhi", "Muon Phi",50,-3.1416,3.1416);
    muonVertex_       = muonDir.make<TH2F>("MuonVertex", "Muon vertex position",50,0,1,100,0,50);

    TFileDirectory genMuonDir = fs_->mkdir("GEN_Muon");
    genMuonMultiplicity_ = genMuonDir.make<TH1F>("genMuonMultiplicity", "Muon Multiplicity", 10,0,10);
    genMuonEnergy_       = genMuonDir.make<TH1F>("genMuonEnergy", "Muon Energy", 100,0,100);
    genMuonPt_           = genMuonDir.make<TH1F>("genMuonPt", "Muon Pt", 100,0,100);
    genMuonEt_           = genMuonDir.make<TH1F>("genMuonEt", "Muon Et", 100,0,100);
    genMuonEta_          = genMuonDir.make<TH1F>("genMuonEta", "Muon Eta",60,-3,3);
    genMuonPhi_          = genMuonDir.make<TH1F>("genMuonPhi", "Muon Phi",50,-3.1416,3.1416);
    genMuonVertex_       = genMuonDir.make<TH2F>("genMuonVertex", "Muon vertex position",50,0,1,100,0,50);

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

    muonEta_->GetXaxis()->SetTitle("eta");
    muonEta_->GetYaxis()->SetTitle("Event Number");
    
    muonPhi_->GetXaxis()->SetTitle("phi");
    muonPhi_->GetYaxis()->SetTitle("Event Number");
    
    muonVertex_->GetXaxis()->SetTitle("V_{z}");
    muonVertex_->GetYaxis()->SetTitle("V_{xy}");

    //------------------------------------------
    genMuonMultiplicity_->GetXaxis()->SetTitle("Number of Muons");
    genMuonMultiplicity_->GetYaxis()->SetTitle("Event Number");
    
    genMuonEnergy_->GetXaxis()->SetTitle("Energy [GeV]");
    genMuonEnergy_->GetYaxis()->SetTitle("Event Number");

    genMuonPt_->GetXaxis()->SetTitle("pT [GeV]");
    genMuonPt_->GetYaxis()->SetTitle("Event Number");

    genMuonEt_->GetXaxis()->SetTitle("ET [GeV]");
    genMuonEt_->GetYaxis()->SetTitle("Event Number");

    genMuonEta_->GetXaxis()->SetTitle("eta");
    genMuonEta_->GetYaxis()->SetTitle("Event Number");

    genMuonPhi_->GetXaxis()->SetTitle("phi");
    genMuonPhi_->GetYaxis()->SetTitle("Event Number");

    genMuonVertex_->GetXaxis()->SetTitle("V_{z}");
    genMuonVertex_->GetYaxis()->SetTitle("V_{xy}");


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
