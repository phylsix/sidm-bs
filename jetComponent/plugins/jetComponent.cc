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
#include "DataFormats/Common/interface/Handle.h"

sidm::jetComponent::jetComponent(const edm::ParameterSet& iConfig):
    patJetToken_(consumes<edm::View<pat::Jet> >(iConfig.getUntrackedParameter<edm::InputTag>("PatJetTag",edm::InputTag("slimmedJets"))) )
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


//
// member functions
//

// ------------ method called for each event  ------------
void
sidm::jetComponent::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
    using namespace edm;
    using namespace std;

    Handle<View<pat::Jet> > jetHdl_;
    iEvent.getByToken(patJetToken_, jetHdl_);
    //const vector<Ptr<pat::Jet> > patJetPtrVec_ = jetHdl_->ptrs();

    /// Loop over jet collection
    int jetNum = 0;
    for (auto iJet = jetHdl_->begin();
              iJet!= jetHdl_->end();
            ++iJet)
    {
        ++jetNum;
    }
    cout<<"Event"<<std::setw(4)<<eventNum_<<":  NumOfJet: "<<jetNum<<endl;
    ++eventNum_;

}


// ------------ method called once each job just before starting event loop  ------------
void 
sidm::jetComponent::beginJob()
{
    patJetTree_ = fs_->make<TTree>("patJetTree", "pat::Jet info");
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
