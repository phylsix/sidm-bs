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

#include "sidm-bs/TreeMaker/interface/TreeMaker.h"

#include <algorithm>

// user include files
#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/Common/interface/Ptr.h"
#include "DataFormats/Math/interface/LorentzVector.h"
#include "DataFormats/Math/interface/deltaR.h"

//
// constants, enums and typedefs
//

//
// static data member definitions
//

//
// constructors and destructor
//
sidm::TreeMaker::TreeMaker(const edm::ParameterSet& iConfig):
        genParticleToken_(consumes<edm::View<reco::GenParticle> >(iConfig.getUntrackedParameter<edm::InputTag>("GenParticleTag",edm::InputTag("prunedGenParticles")))),
        patElectronToken_(consumes<edm::View<pat::Electron> >(iConfig.getUntrackedParameter<edm::InputTag>("PatElectronTag",edm::InputTag("slimmedElectrons")))),
        electronPtLow_(iConfig.getUntrackedParameter<double>("ElectronPtLow")),
        leadElectronPtLow_(iConfig.getUntrackedParameter<double>("LeadElectronPtLow")),
        subleadElectronPtLow_(iConfig.getUntrackedParameter<double>("SubleadElectronPtLow")),
        zMassWindow_(iConfig.getUntrackedParameter<double>("ZMassWindow")),
        realData_(iConfig.getUntrackedParameter<bool>("RealData"))
{
   //now do what ever initialization is needed
   usesResource("TFileService");
   eventNum_ = 0;
}


sidm::TreeMaker::~TreeMaker()
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
sidm::TreeMaker::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  using namespace edm;

  if (!realData_)
  {
    Handle<View<reco::GenParticle> > genParticleHdl_;
    iEvent.getByToken(genParticleToken_, genParticleHdl_);

    const std::vector<edm::Ptr<reco::GenParticle> > genParticlePtrVec_ = genParticleHdl_->ptrs();

    std::vector<edm::Ptr<reco::GenParticle> > genElectronPtrVec_{};
    //std::vector<const reco::Candidate*> electronFromGenZPtrVec_{};

    //* Loop through reco::GenParticles ---------------------- */
    for (std::vector<edm::Ptr<reco::GenParticle> >::const_iterator genParticlePtrIter = genParticlePtrVec_.begin();
         genParticlePtrIter != genParticlePtrVec_.end();
         ++genParticlePtrIter)
    {
        if (abs((*genParticlePtrIter)->pdgId()) != 11 &&
                (*genParticlePtrIter)->pdgId()  != 23) {continue;}

        //* Gen Electron ------------------------------ */
        if (abs((*genParticlePtrIter)->pdgId()) == 11)
        {
            if ((*genParticlePtrIter)->mother()->pdgId() != 23) {continue;}
            if (abs((*genParticlePtrIter)->eta()) > 2.5) {continue;}
            if ((*genParticlePtrIter)->pt() < electronPtLow_) {continue;}
            genElectron_._eventId = eventNum_;
            genElectron_._pt     = (*genParticlePtrIter)->pt();
            genElectron_._eta    = (*genParticlePtrIter)->eta();
            genElectron_._phi    = (*genParticlePtrIter)->phi();
            genElectron_._energy = (*genParticlePtrIter)->energy();

            genElectronTree_->Fill();
            genElectronPtrVec_.push_back(*genParticlePtrIter);
            // std::cout<<"Event# "<<std::setw(5)<<eventNum_<<"    ElectronMother-> "<<(*genParticlePtrIter)->mother()->pdgId()<<std::endl;
        }

        /* Gen Z ----------------------------------
        if ((*genParticlePtrIter)->pdgId() == 23)
        {
            if ((*genParticlePtrIter)->numberOfDaughters() != 2) {continue;}
            const reco::Candidate* d1 = (*genParticlePtrIter)->daughter(0);
            const reco::Candidate* d2 = (*genParticlePtrIter)->daughter(1);
            if ( (abs(d1->pdgId()) != 11) &&
                 (abs(d2->pdgId()) != 11) &&
                 (d1->charge() * d2->charge() != -1) ) {continue;}
            Z_._eventId = eventNum_;
            Z_._mass = (*genParticlePtrIter)->mass();

            ZTree_->Fill();
            electronFromGenZPtrVec_.push_back(d1);
            electronFromGenZPtrVec_.push_back(d2);

            //std::cout<<"Event# "<<std::setw(5)<<eventNum_<<"    Z Mass-->"<<(*genParticlePtrIter)->mass()<<std::endl;
        }
        -----------------------------------------*/
    }

    //* More than 4 Gen Electron ------------------------------ */
    if (genElectronPtrVec_.size() >= 4)
    {
        std::sort(genElectronPtrVec_.begin(), genElectronPtrVec_.end(),
                     [](edm::Ptr<reco::GenParticle>& lhs, edm::Ptr<reco::GenParticle>& rhs)
                     {return lhs->energy() > rhs->energy();} );

        std::vector<std::pair<edm::Ptr<reco::GenParticle>, edm::Ptr<reco::GenParticle> > > ePairInZMassVec_{};
        for (std::vector<edm::Ptr<reco::GenParticle> >::iterator genElectronPtrIter = genElectronPtrVec_.begin();
             genElectronPtrIter != genElectronPtrVec_.end();
             ++genElectronPtrIter)
        {
            if ((*genElectronPtrIter)->pt() < subleadElectronPtLow_) {continue;} // both pt need to be greater than 10GeV
            math::XYZTLorentzVector iP4_= (*genElectronPtrIter)->p4();
            for (std::vector<edm::Ptr<reco::GenParticle> >::iterator jIter = genElectronPtrIter+1;
                 jIter != genElectronPtrVec_.end();
                 ++jIter)
            {
                math::XYZTLorentzVector jP4_ = (*jIter)->p4();
                if ((*genElectronPtrIter)->charge() * (*jIter)->charge() != -1) {continue;} // need to have opposite charge
                if ((iP4_+jP4_).M()<(91.-zMassWindow_) || (iP4_+jP4_).M()>(91+zMassWindow_)) {continue;} // invM need to be within Z +/-10GeV window
                if (std::max((*genElectronPtrIter)->pt(), (*jIter)->pt()) < leadElectronPtLow_) {continue;}
                ePairInZMassVec_.push_back(std::make_pair(*genElectronPtrIter, *jIter));
            }
        }

        //* More than 2 Z reconstructed from Gen Electrons --------------- */
        if(ePairInZMassVec_.size() >= 2)
        {
            std::sort(ePairInZMassVec_.begin(), ePairInZMassVec_.end(),
                         [](std::pair<edm::Ptr<reco::GenParticle>, edm::Ptr<reco::GenParticle> >& lhs,
                            std::pair<edm::Ptr<reco::GenParticle>, edm::Ptr<reco::GenParticle> >& rhs)
                         {
                             float dR_l = deltaR(lhs.first->eta(), lhs.first->phi(), lhs.second->eta(), lhs.second->phi());
                             float dR_r = deltaR(rhs.first->eta(), rhs.first->phi(), rhs.second->eta(), rhs.second->phi());
                             return dR_l < dR_r;
                         });
            genZ1_._eventId = eventNum_;
            genZ1_._mass = (ePairInZMassVec_[0].first->p4() + ePairInZMassVec_[0].second->p4()).M();
            genZ1_._pt   = (ePairInZMassVec_[0].first->p4() + ePairInZMassVec_[0].second->p4()).Pt();

            genZ2_._eventId = eventNum_;
            genZ2_._mass = (ePairInZMassVec_[1].first->p4() + ePairInZMassVec_[1].second->p4()).M();
            genZ2_._pt   = (ePairInZMassVec_[1].first->p4() + ePairInZMassVec_[1].second->p4()).Pt();

            genElectron1_._pt = std::max(ePairInZMassVec_[0].first->pt(), ePairInZMassVec_[0].second->pt());
            genElectron2_._pt = std::min(ePairInZMassVec_[0].first->pt(), ePairInZMassVec_[0].second->pt());
            genElectron3_._pt = std::max(ePairInZMassVec_[1].first->pt(), ePairInZMassVec_[1].second->pt());
            genElectron4_._pt = std::min(ePairInZMassVec_[1].first->pt(), ePairInZMassVec_[1].second->pt());

            genZeeTree_->Fill();

            /*
            std::vector<std::pair<edm::Ptr<reco::GenParticle>, edm::Ptr<reco::GenParticle> > >::iterator electronPairIter;
            for (electronPairIter = ePairInZMassVec_.begin();
                 electronPairIter != ePairInZMassVec_.end();
                 ++electronPairIter)
            {
                float invM_ = (electronPairIter->first->p4() + electronPairIter->second->p4()).M();
                std::cout<<invM_<<" ";
            }
            std::cout<<std::endl;*/
        }
    }
  } // endif (!realData_)

  Handle<View<pat::Electron> > patElectronHdl_;
  iEvent.getByToken(patElectronToken_, patElectronHdl_);

  std::vector<edm::Ptr<pat::Electron> > patElectronPtrVec_ = patElectronHdl_->ptrs();

  //* Got be larger than 4 ------------------------------*/
  if (patElectronPtrVec_.size() >= 4)
  {
      std::sort(patElectronPtrVec_.begin(), patElectronPtrVec_.end(),
                   [](edm::Ptr<pat::Electron>& lhs, edm::Ptr<pat::Electron>& rhs)
                   {return lhs->energy() > rhs->energy();} );
      patElectronPtrVec_.erase(
              std::remove_if(patElectronPtrVec_.begin(), patElectronPtrVec_.end(),
                               [this](const edm::Ptr<pat::Electron>& e)
                               {return (e->pt() < electronPtLow_) || (abs(e->eta()) > 2.5);}),
              patElectronPtrVec_.end()); // In barrel

      //* Loop through pat::Electrons --------------------- */
      std::vector<std::pair<edm::Ptr<pat::Electron>, edm::Ptr<pat::Electron> > > patEPairInZMassVec_{};
      for (std::vector<edm::Ptr<pat::Electron> >::const_iterator patElectronPtrIter = patElectronPtrVec_.begin();
           patElectronPtrIter != patElectronPtrVec_.end();
           ++patElectronPtrIter)
      {
          patElectron_._eventId = eventNum_;
          patElectron_._pt     = (*patElectronPtrIter)->pt();
          patElectron_._eta    = (*patElectronPtrIter)->eta();
          patElectron_._energy = (*patElectronPtrIter)->energy();
          patElectron_._phi    = (*patElectronPtrIter)->phi();

          patElectronTree_->Fill();
          //* --------------------------------------------- */
          if ((*patElectronPtrIter)->pt() < subleadElectronPtLow_) {continue;}
          math::XYZTLorentzVector iP4_{(*patElectronPtrIter)->p4()};
          for (std::vector<edm::Ptr<pat::Electron> >::const_iterator jIter = patElectronPtrIter+1;
               jIter != patElectronPtrVec_.end();
               ++jIter)
          {
              math::XYZTLorentzVector jP4_{(*jIter)->p4()};
              if ((*patElectronPtrIter)->charge() * (*jIter)->charge() != -1) {continue;} // need to have opposite charge
              if ((iP4_+jP4_).M()<(91-zMassWindow_) || (iP4_+jP4_).M()>(91+zMassWindow_)) {continue;} // same as above
              if (std::max((*patElectronPtrIter)->pt(), (*jIter)->pt()) < leadElectronPtLow_) {continue;}
              patEPairInZMassVec_.push_back(std::make_pair(*patElectronPtrIter, *jIter));
          }
      }

      //* More than 2 Z reconstructed from pat::Electrons ---------------*/
      if (patEPairInZMassVec_.size() >= 2)
      {
          std::sort(patEPairInZMassVec_.begin(), patEPairInZMassVec_.end(),
                       [](std::pair<edm::Ptr<pat::Electron>, edm::Ptr<pat::Electron> >& lhs,
                          std::pair<edm::Ptr<pat::Electron>, edm::Ptr<pat::Electron> >& rhs)
                       {
                           float dR_l = deltaR(lhs.first->eta(), lhs.first->phi(), lhs.second->eta(), lhs.second->phi());
                           float dR_r = deltaR(rhs.first->eta(), rhs.first->phi(), rhs.second->eta(), rhs.second->phi());
                           return dR_l < dR_r;
                       });
          patZ1_._eventId = eventNum_;
          patZ1_._mass = (patEPairInZMassVec_[0].first->p4() + patEPairInZMassVec_[0].second->p4()).M();
          patZ1_._pt   = (patEPairInZMassVec_[0].first->p4() + patEPairInZMassVec_[0].second->p4()).Pt();

          patZ2_._eventId = eventNum_;
          patZ2_._mass = (patEPairInZMassVec_[1].first->p4() + patEPairInZMassVec_[1].second->p4()).M();
          patZ2_._pt   = (patEPairInZMassVec_[1].first->p4() + patEPairInZMassVec_[1].second->p4()).Pt();

          patElectron1_._pt = std::max(patEPairInZMassVec_[0].first->pt(), patEPairInZMassVec_[0].second->pt());
          patElectron2_._pt = std::min(patEPairInZMassVec_[0].first->pt(), patEPairInZMassVec_[0].second->pt());
          patElectron3_._pt = std::max(patEPairInZMassVec_[1].first->pt(), patEPairInZMassVec_[1].second->pt());
          patElectron4_._pt = std::min(patEPairInZMassVec_[1].first->pt(), patEPairInZMassVec_[1].second->pt());

          patZeeTree_->Fill();
      }
  }
  eventNum_++;
}


// ------------ method called once each job just before starting event loop  ------------
void 
sidm::TreeMaker::beginJob()
{
  if (!realData_)
  {
    genElectronTree_ = fs_->make<TTree>("genElectronTree", "Gen electrons");
    genElectronTree_->Branch("eventId", &genElectron_._eventId, "eventId/I");
    genElectronTree_->Branch("pt", &genElectron_._pt, "pt/F");
    genElectronTree_->Branch("eta", &genElectron_._eta, "eta/F");
    genElectronTree_->Branch("phi", &genElectron_._phi, "phi/F");
    genElectronTree_->Branch("energy", &genElectron_._energy, "energy/F");

    // ZTree_ = fs_->make<TTree>("ZTree", "Gen Z bosons");
    // ZTree_->Branch("eventId", &Z_._eventId, "eventId/I");
    // ZTree_->Branch("mass", &Z_._mass, "mass/F");

    genZeeTree_ = fs_->make<TTree>("genZeeTree", "Gen info on Zee");
    genZeeTree_->Branch("eventId", &genZ1_._eventId , "eventId/I");
    genZeeTree_->Branch("Z1_Pt", &genZ1_._pt , "Z1_Pt/F");
    genZeeTree_->Branch("Z2_Pt", &genZ2_._pt , "Z2_Pt/F");
    genZeeTree_->Branch("m1ee", &genZ1_._mass , "m1ee/F");
    genZeeTree_->Branch("m2ee", &genZ2_._mass , "m2ee/F");
    genZeeTree_->Branch("e1_Pt", &genElectron1_._pt , "e1_Pt/F");
    genZeeTree_->Branch("e2_Pt", &genElectron2_._pt , "e2_Pt/F");
    genZeeTree_->Branch("e3_Pt", &genElectron3_._pt , "e3_Pt/F");
    genZeeTree_->Branch("e4_Pt", &genElectron4_._pt , "e4_Pt/F");
  }

  patElectronTree_ = fs_->make<TTree>("patElectronTree", "pat::Electron info");
  patElectronTree_->Branch("eventId", &patElectron_._eventId, "eventId/I");
  patElectronTree_->Branch("pt", &patElectron_._pt, "pt/F");
  patElectronTree_->Branch("eta", &patElectron_._eta, "eta/F");
  patElectronTree_->Branch("phi", &patElectron_._phi, "phi/F");
  patElectronTree_->Branch("energy", &patElectron_._energy, "energy/F");

  patZeeTree_ = fs_->make<TTree>("patZeeTree", "pat::Electron info on Zee");
  patZeeTree_->Branch("eventId", &patZ1_._eventId , "eventId/I");
  patZeeTree_->Branch("Z1_Pt", &patZ1_._pt , "Z1_Pt/F");
  patZeeTree_->Branch("Z2_Pt", &patZ2_._pt , "Z2_Pt/F");
  patZeeTree_->Branch("m1ee", &patZ1_._mass , "m1ee/F");
  patZeeTree_->Branch("m2ee", &patZ2_._mass , "m2ee/F");
  patZeeTree_->Branch("e1_Pt", &patElectron1_._pt , "e1_Pt/F");
  patZeeTree_->Branch("e2_Pt", &patElectron2_._pt , "e2_Pt/F");
  patZeeTree_->Branch("e3_Pt", &patElectron3_._pt , "e3_Pt/F");
  patZeeTree_->Branch("e4_Pt", &patElectron4_._pt , "e4_Pt/F");
}

// ------------ method called once each job just after ending the event loop  ------------
void 
sidm::TreeMaker::endJob() 
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
sidm::TreeMaker::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(sidm::TreeMaker);
