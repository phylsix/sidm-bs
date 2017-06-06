// -*- C++ -*-
//
// Package:    sidm-bs/TreeMaker
// Class:      TreeMaker
// 
/**\class TreeMaker TreeMaker.cc sidm-bs/TreeMaker/plugins/TreeMaker.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Weinan Si
//         Created:  Fri, 02 Jun 2017 10:38:41 GMT
//
//


// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/Common/interface/Ptr.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"

#include "TTree.h"

namespace sidm {
    struct Electron{
        int _eventId;
        float _pt;
        float _eta;
        float _phi;
        float _energy;
    };
}
//
// class declaration
//

// If the analyzer does not use TFileService, please remove
// the template argument to the base class so the class inherits
// from  edm::one::EDAnalyzer<> and also remove the line from
// constructor "usesResource("TFileService");"
// This will improve performance in multithreaded jobs.

class TreeMaker : public edm::one::EDAnalyzer<edm::one::SharedResources>  {
// class TreeMaker : public edm::one::EDAnalyzer<>  {
   public:
      explicit TreeMaker(const edm::ParameterSet&);
      ~TreeMaker();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);


   private:
      virtual void beginJob() override;
      virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
      virtual void endJob() override;

      // ----------member data ---------------------------
      edm::Service<TFileService> fs_;
      edm::EDGetTokenT<edm::View<reco::GenParticle> > genParticleToken_;
      TTree* electronTree_;
      sidm::Electron electron_;
      Int_t eventNum_;
};

//
// constants, enums and typedefs
//

//
// static data member definitions
//

//
// constructors and destructor
//
TreeMaker::TreeMaker(const edm::ParameterSet& iConfig):
        genParticleToken_(consumes<edm::View<reco::GenParticle> >(iConfig.getUntrackedParameter<edm::InputTag>("GenParticleTag",edm::InputTag( "prunedGenParticles"))))
{
   //now do what ever initialization is needed
   usesResource("TFileService");
   eventNum_ = 0;
}


TreeMaker::~TreeMaker()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)
    eventNum_ = 0;
}


//
// member functions
//

// ------------ method called for each event  ------------
void
TreeMaker::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   using namespace edm;

   Handle<View<reco::GenParticle> > genParticleHdl_;
   iEvent.getByToken(genParticleToken_, genParticleHdl_);
   const std::vector<edm::Ptr<reco::GenParticle> > genParticlePtrVec_ = genParticleHdl_->ptrs();

   for (std::vector<edm::Ptr<reco::GenParticle> >::const_iterator genParticlePtrIter = genParticlePtrVec_.begin();
        genParticlePtrIter != genParticlePtrVec_.end();
        ++genParticlePtrIter)
   {
       if (abs((*genParticlePtrIter)->pdgId()) != 11) {continue;}
       electron_._eventId = eventNum_;
       electron_._pt = (*genParticlePtrIter)->pt();
       electron_._eta = (*genParticlePtrIter)->eta();
       electron_._phi = (*genParticlePtrIter)->phi();
       electron_._energy = (*genParticlePtrIter)->energy();
       
       electronTree_->Fill();
   }
   eventNum_++;
/* ---------------------
#ifdef THIS_IS_AN_EVENT_EXAMPLE
   Handle<ExampleData> pIn;
   iEvent.getByLabel("example",pIn);
#endif
   
#ifdef THIS_IS_AN_EVENTSETUP_EXAMPLE
   ESHandle<SetupData> pSetup;
   iSetup.get<SetupRecord>().get(pSetup);
#endif
--------------------- */
}


// ------------ method called once each job just before starting event loop  ------------
void 
TreeMaker::beginJob()
{
    electronTree_ = fs_->make<TTree>("ElectronTree", "basic electrons");
    electronTree_->Branch("eventId", &electron_._eventId, "eventId/I");
    electronTree_->Branch("pt", &electron_._pt, "pt/F");
    electronTree_->Branch("eta", &electron_._eta, "eta/F");
    electronTree_->Branch("phi", &electron_._phi, "phi/F");
    electronTree_->Branch("energy", &electron_._energy, "energy/F");
}

// ------------ method called once each job just after ending the event loop  ------------
void 
TreeMaker::endJob() 
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
TreeMaker::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(TreeMaker);
