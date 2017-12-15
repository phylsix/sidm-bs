#ifndef SIDMBS_MUONANALYSIS_H
#define SIDMBS_MUONANALYSIS_H


#include <memory>
#include <vector>

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"
#include "DataFormats/PatCandidates/interface/PackedGenParticle.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/Jet.h"

#include "TTree.h"
#include "TH1F.h"
#include "TH2F.h"
#include "sidm-bs/jetComponent/interface/physicsObject.h"

namespace sidm {

    class MuonAnalysis : public edm::one::EDAnalyzer<edm::one::SharedResources>  {
        public:
            explicit MuonAnalysis(const edm::ParameterSet&);
            ~MuonAnalysis();

            static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);


        private:
            virtual void beginJob() override;
            virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
            virtual void endJob() override;

            // ----------member data ---------------------------

            edm::Service<TFileService> fs_;
            const edm::EDGetTokenT<edm::View<reco::GenParticle> > genParticleTk_;
            const edm::EDGetTokenT<edm::View<pat::PackedGenParticle> > pkdGenTk_;
            const edm::EDGetTokenT<edm::View<pat::PackedCandidate> > pfTk_;
            const edm::EDGetTokenT<edm::View<pat::Electron> > patElectronTk_;
            const edm::EDGetTokenT<edm::View<pat::Muon> > patMuonTk_;
            const edm::EDGetTokenT<edm::View<pat::Jet> > patJetTk_;
            
            TTree* eventTree_;

            TH1F* muonMultiplicity_;
            TH1F* muonEnergy_;
            TH1F* muonPt_;
            TH1F* muonEt_;
            TH1F* muonEta_;
            TH1F* muonPhi_;
            TH2F* muonVertex_;

            TH1F* genMuonMultiplicity_;
            TH1F* genMuonEnergy_;
            TH1F* genMuonPt_;
            TH1F* genMuonEt_;
            TH1F* genMuonEta_;
            TH1F* genMuonPhi_;
            TH2F* genMuonVertex_;
            
            TH1F* matchedMuonMultiplicity_;
            TH1F* matchedMuonEnergy_;
            TH1F* matchedMuonPt_;
            TH1F* matchedMuonEt_;
            TH1F* matchedMuonEta_;
            TH1F* matchedMuonPhi_;
            TH2F* matchedMuonVertex_;
            TH1F* matchedMuonDR_;

            TH1F* unmatchedMuonMultiplicity_;
            TH1F* unmatchedMuonEnergy_;
            TH1F* unmatchedMuonPt_;
            TH1F* unmatchedMuonEt_;
            TH1F* unmatchedMuonEta_;
            TH1F* unmatchedMuonPhi_;
            TH2F* unmatchedMuonVertex_;

            TH1F* matchedDarkphotonMultiplicity_MuMu_;
            TH1F* matchedDarkphotonInvm_MuMu_;
            TH1F* matchedDarkphotonDR_MuMu_;
            Int_t eventNum_;
    };

}  // namespace sidm



#endif
