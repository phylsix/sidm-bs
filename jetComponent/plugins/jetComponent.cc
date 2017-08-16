// -*- C++ -*-
//
// Package:    sidm-bs/jetComponent
// Class:      jetComponent
// 
/**\class jetComponent jetComponent.cc sidm-bs/jetComponent/plugins/jetComponent.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  weinan si
//         Created:  Thu, 10 Aug 2017 23:23:35 GMT
//
//

#include "sidm-bs/jetComponent/interface/jetComponent.h"
#include "sidm-bs/jetComponent/interface/physicsObject.h"
#include "DataFormats/Common/interface/Handle.h"

sidm::jetComponent::jetComponent(const edm::ParameterSet& iConfig):
    genParticleTk_(consumes<edm::View<reco::GenParticle> >(iConfig.getUntrackedParameter<edm::InputTag>("GenParticleTag_", edm::InputTag("prunedGenParticles")))),
    genJetTk_(consumes<edm::View<reco::GenJet> >(iConfig.getUntrackedParameter<edm::InputTag>("GenJetTag_1", edm::InputTag("slimmedGenJets")))),
    genJetAK8Tk_(consumes<edm::View<reco::GenJet> >(iConfig.getUntrackedParameter<edm::InputTag>("GenJetTag_2", edm::InputTag("slimmedGenJetsAK8")))),
    patJetTk_(consumes<edm::View<pat::Jet> >(iConfig.getUntrackedParameter<edm::InputTag>("PatJetTag_1",edm::InputTag("slimmedJets"))) ),
    patJetAK8Tk_(consumes<edm::View<pat::Jet> >(iConfig.getUntrackedParameter<edm::InputTag>("PatJetTag_2",edm::InputTag("slimmedJetsAK8"))) ),
    patJetPuppiTk_(consumes<edm::View<pat::Jet> >(iConfig.getUntrackedParameter<edm::InputTag>("PatJetTag_3",edm::InputTag("slimmedJetsPuppi"))) ),
    patJetAK8CHSTk_(consumes<edm::View<pat::Jet> >(iConfig.getUntrackedParameter<edm::InputTag>("PatJetTag_4",edm::InputTag("slimmedJetsAK8PFCHSSoftDropPacked"))) ),
    patJetAK8PuppiTk_(consumes<edm::View<pat::Jet> >(iConfig.getUntrackedParameter<edm::InputTag>("PatJetTag_5",edm::InputTag("slimmedJetsAK8PFPuppiSoftDropPacked"))) )
{
    //now do what ever initialization is needed
    usesResource("TFileService");
    eventNum_ = 0;
}


sidm::jetComponent::~jetComponent()
{

    // do anything here that needs to be done at desctruction time
    // (e.g. close files, deallocate resources etc.)
    eventNum_ = 0;

}



// ------------ method called for each event  ------------
void
sidm::jetComponent::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
    using namespace edm;
    using namespace std;

    Handle<View<reco::GenJet> > genJetHdl_;
    iEvent.getByToken(genJetTk_, genJetHdl_);

    for ( const auto& j : *genJetHdl_ ) {
        jetGen_._eventId = eventNum_;
        jetGen_._eta = j.eta();
        jetGen_._phi = j.phi();
        jetGen_._pt = j.pt();
        jetGen_._energy = j.energy();
        jetGen_._mass = j.mass();
        Gen_slimmedGenJets_->Fill();
    }

    Handle<View<reco::GenParticle> > genParticleHdl_;
    iEvent.getByToken(genParticleTk_, genParticleHdl_);

    for ( const auto& p : *genParticleHdl_ ) {
        if ( abs(p.pdgId()) != 11 ) {continue;}
        epGen_._eventId = eventNum_;
        epGen_._eta = p.eta();
        epGen_._phi = p.phi();
        epGen_._pt = p.pt();
        epGen_._energy = p.energy();
        Gen_electrons_->Fill();
    }


    Handle<View<reco::GenJet> > genJetAK8Hdl_;
    iEvent.getByToken(genJetAK8Tk_, genJetAK8Hdl_);

    for ( const auto& j : *genJetAK8Hdl_ ) {
        jetGenAK8_._eventId = eventNum_;
        jetGenAK8_._eta = j.eta();
        jetGenAK8_._phi = j.phi();
        jetGenAK8_._pt = j.pt();
        jetGenAK8_._energy = j.energy();
        jetGenAK8_._mass = j.mass();
        Gen_slimmedGenJetsAK8_->Fill();
    }


    Handle<View<pat::Jet> > patJetHdl_;
    iEvent.getByToken(patJetTk_, patJetHdl_);
    //const vector<Ptr<pat::Jet> > patJetPtrVec_ = jetHdl_->ptrs();

    for ( const auto& j : *patJetHdl_ ) {
        jetPat_._eventId = eventNum_;
        jetPat_._eta = j.eta();
        jetPat_._phi = j.phi();
        jetPat_._pt = j.pt();
        jetPat_._energy = j.energy();
        jetPat_._mass = j.mass();
        Pat_slimmedJets_->Fill();
    }

    Handle<View<pat::Jet> > patJetAK8Hdl_;
    iEvent.getByToken(patJetAK8Tk_, patJetAK8Hdl_);

    for ( const auto& j: *patJetAK8Hdl_ ) {
        jetPatAK8_._eventId = eventNum_;
        jetPatAK8_._eta = j.eta();
        jetPatAK8_._phi = j.phi();
        jetPatAK8_._pt = j.pt();
        jetPatAK8_._energy = j.energy();
        jetPatAK8_._mass = j.mass();
        Pat_slimmedJetsAK8_->Fill();
    }

    Handle<View<pat::Jet> > patJetPuppiHdl_;
    iEvent.getByToken(patJetPuppiTk_, patJetPuppiHdl_);

    for ( const auto& j: *patJetPuppiHdl_ ) {
        jetPatPuppi_._eventId = eventNum_;
        jetPatPuppi_._eta = j.eta();
        jetPatPuppi_._phi = j.phi();
        jetPatPuppi_._pt = j.pt();
        jetPatPuppi_._energy = j.energy();
        jetPatPuppi_._mass = j.mass();
        Pat_slimmedJetsPuppi_->Fill();
    }

    Handle<View<pat::Jet> > patJetAK8CHSHdl_;
    iEvent.getByToken(patJetAK8CHSTk_, patJetAK8CHSHdl_);

    for ( const auto& j: *patJetAK8CHSHdl_ ) {
        jetPatAK8PFCHSSoftDropPacked_._eventId = eventNum_;
        jetPatAK8PFCHSSoftDropPacked_._eta = j.eta();
        jetPatAK8PFCHSSoftDropPacked_._phi = j.phi();
        jetPatAK8PFCHSSoftDropPacked_._pt = j.pt();
        jetPatAK8PFCHSSoftDropPacked_._energy = j.energy();
        jetPatAK8PFCHSSoftDropPacked_._mass = j.mass();
        Pat_slimmedJetsAK8PFCHSSoftDropPacked_->Fill();
    }

    Handle<View<pat::Jet> > patJetAK8PuppiHdl_;
    iEvent.getByToken(patJetAK8PuppiTk_, patJetAK8PuppiHdl_);

    for ( const auto& j: *patJetAK8PuppiHdl_ ) {
        jetPatAK8PFPuppiSoftDropPacked_._eventId = eventNum_;
        jetPatAK8PFPuppiSoftDropPacked_._eta = j.eta();
        jetPatAK8PFPuppiSoftDropPacked_._phi = j.phi();
        jetPatAK8PFPuppiSoftDropPacked_._pt = j.pt();
        jetPatAK8PFPuppiSoftDropPacked_._energy = j.energy();
        jetPatAK8PFPuppiSoftDropPacked_._mass = j.mass();
        Pat_slimmedJetsAK8PFPuppiSoftDropPacked_->Fill();
    }


    /// Loop over jet collection
    //int jetNum = 0;
    //for (auto iJet = jetHdl_->begin();
    //          iJet!= jetHdl_->end();
    //        ++iJet)
    //{
    //    ++jetNum;
    //}
    //cout<<"Event"<<std::setw(4)<<eventNum_<<":  NumOfJet: "<<jetNum<<endl;
    ++eventNum_;

}


// ------------ method called once each job just before starting event loop  ------------
void 
sidm::jetComponent::beginJob()
{
    /// generated electrons
    Gen_electrons_ = fs_->make<TTree>("Gen_electrons", "gen electrons");

    Gen_electrons_->Branch("eventId", &epGen_._eventId, "eventId/I");
    Gen_electrons_->Branch("eta",     &epGen_._eta,     "eta/F");
    Gen_electrons_->Branch("phi",     &epGen_._phi,     "phi/F");
    Gen_electrons_->Branch("pt",      &epGen_._pt,      "pt/F");
    Gen_electrons_->Branch("energy",  &epGen_._energy,  "energy/F");


    /// generated jets
    Gen_slimmedGenJets_ = fs_->make<TTree>("Gen_slimmedGenJets", "slimmedGenJets");

    Gen_slimmedGenJets_->Branch("eventId", &jetGen_._eventId, "eventId/I");
    Gen_slimmedGenJets_->Branch("eta",     &jetGen_._eta,     "eta/F");
    Gen_slimmedGenJets_->Branch("phi",     &jetGen_._phi,     "phi/F");
    Gen_slimmedGenJets_->Branch("pt",      &jetGen_._pt,      "pt/F");
    Gen_slimmedGenJets_->Branch("energy",  &jetGen_._energy,  "energy/F");
    Gen_slimmedGenJets_->Branch("mass",    &jetGen_._mass,    "mass/F");

    Gen_slimmedGenJetsAK8_ = fs_->make<TTree>("Gen_slimmedGenJetsAK8", "slimmedGenJetsAK8");

    Gen_slimmedGenJetsAK8_->Branch("eventId", &jetGenAK8_._eventId, "eventId/I");
    Gen_slimmedGenJetsAK8_->Branch("eta",     &jetGenAK8_._eta,     "eta/F");
    Gen_slimmedGenJetsAK8_->Branch("phi",     &jetGenAK8_._phi,     "phi/F");
    Gen_slimmedGenJetsAK8_->Branch("pt",      &jetGenAK8_._pt,      "pt/F");
    Gen_slimmedGenJetsAK8_->Branch("energy",  &jetGenAK8_._energy,  "energy/F");
    Gen_slimmedGenJetsAK8_->Branch("mass",    &jetGenAK8_._mass,    "mass/F");


    /// reconstructed jets
    Pat_slimmedJets_ = fs_->make<TTree>("PAT_slimmedJets", "slimmedJets");

    Pat_slimmedJets_->Branch("eventId", &jetPat_._eventId, "eventId/I");
    Pat_slimmedJets_->Branch("eta",     &jetPat_._eta,     "eta/F");
    Pat_slimmedJets_->Branch("phi",     &jetPat_._phi,     "phi/F");
    Pat_slimmedJets_->Branch("pt",      &jetPat_._pt,      "pt/F");
    Pat_slimmedJets_->Branch("energy",  &jetPat_._energy,  "energy/F");
    Pat_slimmedJets_->Branch("mass",    &jetPat_._mass,    "mass/F");

    Pat_slimmedJetsAK8_ = fs_->make<TTree>("Pat_slimmedJetsAK8", "slimmedJetsAK8");

    Pat_slimmedJetsAK8_->Branch("eventId", &jetPatAK8_._eventId, "eventId/I");
    Pat_slimmedJetsAK8_->Branch("eta",     &jetPatAK8_._eta,     "eta/F");
    Pat_slimmedJetsAK8_->Branch("phi",     &jetPatAK8_._phi,     "phi/F");
    Pat_slimmedJetsAK8_->Branch("pt",      &jetPatAK8_._pt,      "pt/F");
    Pat_slimmedJetsAK8_->Branch("energy",  &jetPatAK8_._energy,  "energy/F");
    Pat_slimmedJetsAK8_->Branch("mass",    &jetPatAK8_._mass,    "mass/F");

    Pat_slimmedJetsPuppi_ = fs_->make<TTree>("Pat_slimmedJetsPuppi", "slimmedJetsPuppi");

    Pat_slimmedJetsPuppi_->Branch("eventId", &jetPatPuppi_._eventId, "eventId/I");
    Pat_slimmedJetsPuppi_->Branch("eta",     &jetPatPuppi_._eta,     "eta/F");
    Pat_slimmedJetsPuppi_->Branch("phi",     &jetPatPuppi_._phi,     "phi/F");
    Pat_slimmedJetsPuppi_->Branch("pt",      &jetPatPuppi_._pt,      "pt/F");
    Pat_slimmedJetsPuppi_->Branch("energy",  &jetPatPuppi_._energy,  "energy/F");
    Pat_slimmedJetsPuppi_->Branch("mass",    &jetPatPuppi_._mass,    "mass/F");

    Pat_slimmedJetsAK8PFCHSSoftDropPacked_ = fs_->make<TTree>("Pat_slimmedJetsAK8PFCHSSoftDropPacked", "slimmedJetsAK8PFCHSSoftDropPacked");

    Pat_slimmedJetsAK8PFCHSSoftDropPacked_->Branch("eventId", &jetPatAK8PFCHSSoftDropPacked_._eventId, "eventId/I");
    Pat_slimmedJetsAK8PFCHSSoftDropPacked_->Branch("eta",     &jetPatAK8PFCHSSoftDropPacked_._eta,     "eta/F");
    Pat_slimmedJetsAK8PFCHSSoftDropPacked_->Branch("phi",     &jetPatAK8PFCHSSoftDropPacked_._phi,     "phi/F");
    Pat_slimmedJetsAK8PFCHSSoftDropPacked_->Branch("pt",      &jetPatAK8PFCHSSoftDropPacked_._pt,      "pt/F");
    Pat_slimmedJetsAK8PFCHSSoftDropPacked_->Branch("energy",  &jetPatAK8PFCHSSoftDropPacked_._energy,  "energy/F");
    Pat_slimmedJetsAK8PFCHSSoftDropPacked_->Branch("mass",    &jetPatAK8PFCHSSoftDropPacked_._mass,    "mass/F");

    Pat_slimmedJetsAK8PFPuppiSoftDropPacked_ = fs_->make<TTree>("Pat_slimmedJetsAK8PFPuppiSoftDropPacked", "slimmedJetsAK8PFPuppiSoftDropPacked");

    Pat_slimmedJetsAK8PFPuppiSoftDropPacked_->Branch("eventId", &jetPatAK8PFPuppiSoftDropPacked_._eventId, "eventId/I");
    Pat_slimmedJetsAK8PFPuppiSoftDropPacked_->Branch("eta",     &jetPatAK8PFPuppiSoftDropPacked_._eta,     "eta/F");
    Pat_slimmedJetsAK8PFPuppiSoftDropPacked_->Branch("phi",     &jetPatAK8PFPuppiSoftDropPacked_._phi,     "phi/F");
    Pat_slimmedJetsAK8PFPuppiSoftDropPacked_->Branch("pt",      &jetPatAK8PFPuppiSoftDropPacked_._pt,      "pt/F");
    Pat_slimmedJetsAK8PFPuppiSoftDropPacked_->Branch("energy",  &jetPatAK8PFPuppiSoftDropPacked_._energy,  "energy/F");
    Pat_slimmedJetsAK8PFPuppiSoftDropPacked_->Branch("mass",    &jetPatAK8PFPuppiSoftDropPacked_._mass,    "mass/F");

}

// ------------ method called once each job just after ending the event loop  ------------
void 
sidm::jetComponent::endJob() 
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
sidm::jetComponent::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(sidm::jetComponent);
