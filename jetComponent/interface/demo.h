#ifndef SIDMBS_DEMO_H
#define SIDMBS_DEMO_H


#include <memory>
#include <vector>

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/Jet.h"

#include "TTree.h"
#include "sidm-bs/jetComponent/interface/physicsObject.h"

namespace sidm {

    class demo : public edm::one::EDAnalyzer<edm::one::SharedResources>  {
        public:
            explicit demo(const edm::ParameterSet&);
            ~demo();

            static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);


        private:
            virtual void beginJob() override;
            virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
            virtual void endJob() override;

            // ----------member data ---------------------------

            edm::Service<TFileService> fs_;
            const edm::EDGetTokenT<edm::View<reco::GenParticle> > genParticleTk_;
            const edm::EDGetTokenT<edm::View<pat::PackedCandidate> > pfTk_;
            const edm::EDGetTokenT<edm::View<pat::Electron> > patElectronTk_;
            const edm::EDGetTokenT<edm::View<pat::Jet> > patJetTk_;
            
            TTree* eventTree_;


            Int_t eventNum_;
    };

}  // namespace sidm



#endif
