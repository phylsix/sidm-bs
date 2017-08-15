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
#include "DataFormats/JetReco/interface/PFJet.h"
#include "DataFormats/JetReco/interface/PFJetCollection.h"
#include "DataFormats/PatCandidates/interface/PackedGenParticle.h"
#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"

#include "TTree.h"
#include "sidm-bs/jetComponent/interface/physicsObject.h"

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
            const edm::EDGetTokenT<edm::View<reco::GenJet> > genJetTk_;
            const edm::EDGetTokenT<edm::View<reco::GenJet> > genJetAK8Tk_;
            const edm::EDGetTokenT<edm::View<pat::Jet> > patJetTk_;
            const edm::EDGetTokenT<edm::View<pat::Jet> > patJetAK8Tk_;
            const edm::EDGetTokenT<edm::View<pat::Jet> > patJetPuppiTk_;
            const edm::EDGetTokenT<edm::View<pat::Jet> > patJetAK8CHSTk_;
            const edm::EDGetTokenT<edm::View<pat::Jet> > patJetAK8PuppiTk_;
            const edm::EDGetTokenT<edm::View<reco::GenParticle> > genParticleTk_;

            TTree* Gen_slimmedGenJets_;
            TTree* Gen_slimmedGenJetsAK8_;

            sidm::Jet jetGen_;
            sidm::Jet jetGenAK8_;

            TTree* Pat_slimmedJets_;
            TTree* Pat_slimmedJetsAK8_;
            TTree* Pat_slimmedJetsPuppi_;
            TTree* Pat_slimmedJetsAK8PFCHSSoftDropPacked_;
            TTree* Pat_slimmedJetsAK8PFPuppiSoftDropPacked_;

            sidm::Jet jetPat_;
            sidm::Jet jetPatAK8_;
            sidm::Jet jetPatPuppi_;
            sidm::Jet jetPatAK8PFCHSSoftDropPacked_;
            sidm::Jet jetPatAK8PFPuppiSoftDropPacked_;

            Int_t eventNum_;
    };

}

#endif
