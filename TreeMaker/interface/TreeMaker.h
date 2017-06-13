#ifndef SIDMBS_TREEMAKER_TREEMAKER_H
#define SIDMBS_TREEMAKER_TREEMAKER_H

#include <memory>
#include <vector>

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "DataFormats/Common/interface/Ptr.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/PatCandidates/interface/Electron.h"

#include "TTree.h"

namespace sidm {
    struct Electron{
        int _eventId;
        float _pt;
        float _eta;
        float _phi;
        float _energy;
    };
    struct Z{
        int _eventId;
        float _mass;
        float _pt;
    };
    
    class TreeMaker : public edm::one::EDAnalyzer<edm::one::SharedResources>  {
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
          edm::EDGetTokenT<edm::View<pat::Electron> > patElectronToken_;
          double electronPtLow_;
          double leadElectronPtLow_;
          double subleadElectronPtLow_;
          double zMassWindow_;
          bool realData_;
    
          TTree* genElectronTree_;
          sidm::Electron genElectron_;
    
          // TTree* ZTree_;
          // sidm::Z Z_;
    
          TTree* genZeeTree_;
          sidm::Z genZ1_;
          sidm::Z genZ2_;
          sidm::Electron genElectron1_;
          sidm::Electron genElectron2_;
          sidm::Electron genElectron3_;
          sidm::Electron genElectron4_;
    
          TTree* patElectronTree_;
          sidm::Electron patElectron_;
    
          TTree* patZeeTree_;
          sidm::Z patZ1_;
          sidm::Z patZ2_;
          sidm::Electron patElectron1_;
          sidm::Electron patElectron2_;
          sidm::Electron patElectron3_;
          sidm::Electron patElectron4_;
    
          Int_t eventNum_;
    };

}

#endif  // SIDMBS_TREEMAKER_TREEMAKER_H_
