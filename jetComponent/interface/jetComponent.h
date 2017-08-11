#ifndef SIDMBS_JETCOMPONENT_H
#define SIDMBS_JETCOMPONENT_H

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
#include "DataFormats/PatCandidates/interface/Jet.h"

#include "TTree.h"

namespace sidm {

    class jetComponent : public edm::one::EDAnalyzer<edm::one::SharedResources>  {
        public:
            explicit jetComponent(const edm::ParameterSet&);
            ~jetComponent();

            static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);


        private:
            virtual void beginJob() override;
            virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
            virtual void endJob() override;

            // ----------member data ---------------------------

            edm::Service<TFileService> fs_;
            const edm::EDGetTokenT<edm::View<pat::Jet> > patJetToken_;

            TTree* patJetTree_;

            Int_t eventNum_;
    };

}

#endif
