#ifndef SIDMBS_MCVALIDATION_H
#define SIDMBS_MCVALIDATION_H

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
#include "DataFormats/PatCandidates/interface/PackedGenParticle.h"
#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"

#include "TTree.h"
#include "sidm-bs/jetComponent/interface/physicsObject.h"

namespace sidm {

    class mcValidation : public edm::one::EDAnalyzer<edm::one::SharedResources>  {
        public:
            explicit mcValidation(const edm::ParameterSet&);
            ~mcValidation();

            static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);


        private:
            virtual void beginJob() override;
            virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
            virtual void endJob() override;

            // ----------member data ---------------------------

            edm::Service<TFileService> fs_;
            const edm::EDGetTokenT<edm::View<reco::GenParticle> > genParticleTk_;
            
            TTree* eventTree_;

            TTree* darkPhoton_reco_;
            sidm::Zp zp_;

            TTree* pscalar_reco_;
            sidm::Ps ps_;

            Int_t eventNum_;
            Int_t electron_N;
            Int_t positron_N;
            Int_t electron_from_zp_N;
            Int_t positron_from_zp_N;
            Int_t zp_N;
            Int_t ps_N;
    };

}

#endif
