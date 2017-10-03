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
    pkdGenTk_(consumes<edm::View<pat::PackedGenParticle> >(iConfig.getUntrackedParameter<edm::InputTag>("PatPkdGenTag_", edm::InputTag("packedGenParticles")))),
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

    Handle<View<pat::PackedGenParticle> > pkdGenHdl_;
    iEvent.getByToken(pkdGenTk_, pkdGenHdl_);

    Handle<View<pat::Electron> > patElectronHdl_;
    iEvent.getByToken(patElectronTk_, patElectronHdl_);

    Handle<View<pat::Jet> > patJetHdl_;
    iEvent.getByToken(patJetTk_, patJetHdl_);

    Handle<View<pat::PackedCandidate> > pfHdl_;
    iEvent.getByToken(pfTk_, pfHdl_);

    pkdGenElectron_N = 0; /// number of electrons who are coming from darkphotons
    pfElectron_N     = 0;
    pfGamma_N        = 0;

    for (const auto& zp : *genParticleHdl_) {
        if (zp.pdgId() != 32) continue;
        //vector<const pat::PackedGenParticle*> electronsFromSingelDarkPhoton{};
        vector<const reco::GenParticle*> electronsFromSingelDarkPhoton{};

        //for (const auto& p : *pkdGenHdl_) {
        for (const auto& p : *genParticleHdl_) {
            if (p.status() != 1 || abs(p.pdgId()) != 11) continue;
            const reco::Candidate* motherInPrunedCollection(p.mother(0));
            if (motherInPrunedCollection != nullptr && sidm::is_ancestor(&zp, motherInPrunedCollection)) {
                ++pkdGenElectron_N;
                electronsFromSingelDarkPhoton.push_back(&p);
            }
        }

        if (electronsFromSingelDarkPhoton.size()!=2) {
            cout<<"Event"<<eventNum_<<": Darkphoton is not mother of two electrons, instead- "<<electronsFromSingelDarkPhoton.size()<<endl;
            continue; /// This darkphoton does not have two electron daughters
        }

        if (electronsFromSingelDarkPhoton[0]->charge() > 0 && electronsFromSingelDarkPhoton[1]->charge() < 0) {
            darkPhotonFromGenEle_ = sidm::Zp(electronsFromSingelDarkPhoton[1], electronsFromSingelDarkPhoton[0]);
        } else if (electronsFromSingelDarkPhoton[0]->charge() < 0 && electronsFromSingelDarkPhoton[1]->charge() > 0) {
            darkPhotonFromGenEle_ = sidm::Zp(electronsFromSingelDarkPhoton[0], electronsFromSingelDarkPhoton[1]);
        } else {
            cout<<"Event"<<eventNum_<<": Darkphoton is mother of two electrons, but not opposite charge, skip-";
            continue;
        }

        darkPhotonFromGenEle_._eventId = eventNum_;
        darkPhotonFromGenEle_._pt      = zp.pt();
        darkPhotonFromGenEle_._eta     = zp.eta();
        darkPhotonFromGenEle_._mass    = zp.mass();

        darkPhotonFromGenElectronsTree_->Fill();
        //assert(electronsFromSingelDarkPhoton.size() == 2);
    }
    //pkdGenElectron_N = count_if(cbegin(*pkdGenHdl_), cend(*pkdGenHdl_),
    //        [](const auto& p){ return abs(p.pdgId()) == 11 && p.mother(0)->pdgId() == 32; });
    pfElectron_N = count_if(cbegin(*pfHdl_), cend(*pfHdl_),
            [](const auto& p){ return abs(p.pdgId()) == 11; });
    pfGamma_N = count_if(cbegin(*pfHdl_), cend(*pfHdl_),
            [](const auto& p){ return p.pdgId() == 22; });
    eventTree_->Fill();

    ++eventNum_;
}


// ------------ method called once each job just before starting event loop  ------------
void 
sidm::electronFinder::beginJob()
{
    /// Statistics of multiplicities of collections per event
    eventTree_ = fs_->make<TTree>("eventTree", "information per event");
    
    //eventTree_->Branch("numOfEleInPrunedGen", &genElectron_N, "numOfEleInPrunedGen/I");
    eventTree_->Branch("numOfEleInPackedGen", &pkdGenElectron_N, "numOfEleInPackedGen/I");
    eventTree_->Branch("numOfEleInPfCands", &pfElectron_N, "numOfEleInPfCands/I");
    eventTree_->Branch("numOfPhoInPfCands", &pfGamma_N, "numOfPhoInPfCands/I");
    // eventTree_->Branch("numberOfElectrons", &electron_N, "numberOfElectrons/I");

    darkPhotonFromGenElectronsTree_ = fs_->make<TTree>("darkPhotonFromGenElectrons",
            "darkphotons reco from electrons in packedGenParticles collection");
    darkPhotonFromGenElectronsTree_->Branch("eventId", &darkPhotonFromGenEle_._eventId, "eventId/I");
    darkPhotonFromGenElectronsTree_->Branch("pt",      &darkPhotonFromGenEle_._pt, "pt/F");
    darkPhotonFromGenElectronsTree_->Branch("eta",     &darkPhotonFromGenEle_._eta, "eta/F");
    darkPhotonFromGenElectronsTree_->Branch("mass",    &darkPhotonFromGenEle_._mass, "mass/F");
    darkPhotonFromGenElectronsTree_->Branch("invm",    &darkPhotonFromGenEle_._invM, "invm/F");
    darkPhotonFromGenElectronsTree_->Branch("pt_e",    &darkPhotonFromGenEle_.e._pt, "pt_e/F");
    darkPhotonFromGenElectronsTree_->Branch("pt_p",    &darkPhotonFromGenEle_.p._pt, "pt_p/F");
    darkPhotonFromGenElectronsTree_->Branch("dR_ep",   &darkPhotonFromGenEle_._dR, "dR_ep/F");
    darkPhotonFromGenElectronsTree_->Branch("dEta_ep", &darkPhotonFromGenEle_._dEta, "dEta_ep/F");
    darkPhotonFromGenElectronsTree_->Branch("dPhi_ep", &darkPhotonFromGenEle_._dPhi, "dPhi_ep/F");
    darkPhotonFromGenElectronsTree_->Branch("dv_x",    &darkPhotonFromGenEle_._dv_x, "dv_x/F");
    darkPhotonFromGenElectronsTree_->Branch("dv_y",    &darkPhotonFromGenEle_._dv_y, "dv_y/F");
    darkPhotonFromGenElectronsTree_->Branch("dv_z",    &darkPhotonFromGenEle_._dv_z, "dv_z/F");


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
