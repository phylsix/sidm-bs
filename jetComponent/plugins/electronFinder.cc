// -*- C++ -*-
//
// Package:    sidm-bs/electronFinder
// Class:      electronFinder
#include <algorithm>
#include <cmath>

#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/Common/interface/Ptr.h"
#include "DataFormats/Math/interface/LorentzVector.h"
#include "DataFormats/Math/interface/deltaR.h"

#include "sidm-bs/jetComponent/interface/physicsObject.h"
#include "sidm-bs/jetComponent/interface/utilities.h"

#include "sidm-bs/jetComponent/interface/electronFinder.h"

sidm::electronFinder::electronFinder(const edm::ParameterSet& iConfig):
    genParticleTk_(consumes<edm::View<reco::GenParticle> >(iConfig.getUntrackedParameter<edm::InputTag>("GenParticleTag_", edm::InputTag("prunedGenParticles")))),
    pfTk_(consumes<edm::View<pat::PackedCandidate> >(iConfig.getUntrackedParameter<edm::InputTag>("PatPfTag_", edm::InputTag("packedPFCandidates")))),
    pkdGenTk_(consumes<edm::View<pat::PackedGenParticle> >(iConfig.getUntrackedParameter<edm::InputTag>("PatPkdGenTag_", edm::InputTag("packedGenParticles")))),
    patElectronTk_(consumes<edm::View<pat::Electron> >(iConfig.getUntrackedParameter<edm::InputTag>("PatElectronTag_", edm::InputTag("slimmedElectrons")))),
    patJetTk_(consumes<edm::View<pat::Jet> >(iConfig.getUntrackedParameter<edm::InputTag>("PatJetTag_", edm::InputTag("slimmedJets"))))
{
    usesResource("TFileService");
    eventNum_ = 0;
}


sidm::electronFinder::~electronFinder()
{

    std::cout << "\n\nNumber of events: " << eventNum_ << "\n\n";
    eventNum_ = 0;

}

// ------------ method called for each event  ------------
void
sidm::electronFinder::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
    using namespace edm;
    using namespace std;

    Handle<View<reco::GenParticle> > genParticleHdl_;
    iEvent.getByToken(genParticleTk_, genParticleHdl_);

    Handle<View<pat::PackedGenParticle> > pkdGenHdl_;
    iEvent.getByToken(pkdGenTk_, pkdGenHdl_);

    Handle<View<pat::Electron> > patElectronHdl_;
    iEvent.getByToken(patElectronTk_, patElectronHdl_);

    Handle<View<pat::Jet> > patJetHdl_;
    iEvent.getByToken(patJetTk_, patJetHdl_);

    Handle<View<pat::PackedCandidate> > pfHdl_;
    iEvent.getByToken(pfTk_, pfHdl_);

    pkdGenElectron_N = 0; /// number of electrons who are coming from darkphotons
    pfElectron_N     = 0;
    pfGamma_N        = 0;
    pkdGenEleCutApplied_N = 0;
    epOfZpInJet_N = 0;

    vector<Ptr<pat::PackedGenParticle> > pkdGenPtr_ = pkdGenHdl_->ptrs();
    vector<Ptr<pat::PackedCandidate> >   pfPtr_     = pfHdl_->ptrs();
    vector<Ptr<pat::Jet> >               patJetPtr_ = patJetHdl_->ptrs();
    
    vector<sidm::Zp> darkPhotonFromElectronsInPackedGen{};
    for (const auto& zp : *genParticleHdl_) {
        if (zp.pdgId() != 32) continue;
        vector<Ptr<pat::PackedGenParticle> > electronsFromSingleDarkPhoton{};
        electronsFromSingleDarkPhoton.clear();

        for (const auto& p : pkdGenPtr_) {
            if (p->status() != 1 || abs(p->pdgId()) != 11) continue;
            const reco::Candidate* motherInPrunedCollection(p->mother(0));
            if (motherInPrunedCollection != nullptr && sidm::is_ancestor(&zp, motherInPrunedCollection)) {
                ++pkdGenElectron_N;
                electronsFromSingleDarkPhoton.push_back(p);
            }
        }

        if (electronsFromSingleDarkPhoton.size()!=2) {
            cout<<"Event"<<eventNum_<<": Darkphoton is not mother of two electrons, instead- "<<electronsFromSingleDarkPhoton.size()<<endl;
            continue; /// This darkphoton does not have two electron daughters
        }

        if (electronsFromSingleDarkPhoton[0]->charge() > 0 && electronsFromSingleDarkPhoton[1]->charge() < 0) {
            //sidm::Ep ele(electronsFromSingleDarkPhoton[1]), pos(electronsFromSingleDarkPhoton[0]);
            //darkPhotonFromGenEle_ = sidm::Zp(ele, pos);
            darkPhotonFromGenEle_ = sidm::Zp(electronsFromSingleDarkPhoton[1], electronsFromSingleDarkPhoton[0]);
        } else if (electronsFromSingleDarkPhoton[0]->charge() < 0 && electronsFromSingleDarkPhoton[1]->charge() > 0) {
            //sidm::Ep ele(electronsFromSingleDarkPhoton[0]), pos(electronsFromSingleDarkPhoton[1]);
            //darkPhotonFromGenEle_ = sidm::Zp(ele, pos);
            darkPhotonFromGenEle_ = sidm::Zp(electronsFromSingleDarkPhoton[0], electronsFromSingleDarkPhoton[1]);
        } else {
            cout<<"Event"<<eventNum_<<": Darkphoton is mother of two electrons, but not opposite charge, skip-\n";
            continue;
        }

        darkPhotonFromGenEle_.setVertex(zp.daughter(0)->vertex()); // complete with vertex info which lack from packedGenParticle

        darkPhotonFromGenEle_._eventId = eventNum_;
        darkPhotonFromGenEle_._pt      = zp.pt();
        darkPhotonFromGenEle_._eta     = zp.eta();
        darkPhotonFromGenEle_._mass    = zp.mass();

        darkPhotonFromGenElectronsTree_->Fill();
        darkPhotonFromElectronsInPackedGen.push_back(darkPhotonFromGenEle_);
    }
    // cout<<"Electron in darkPhotonFromElectronsInPackedGen[0]: eta:"<<darkPhotonFromElectronsInPackedGen[0].e._eta<<" phi:"<<darkPhotonFromElectronsInPackedGen[0].e._phi<<endl;
    // cout<<"Electron in darkPhotonFromElectronsInPackedGen[1]: eta:"<<darkPhotonFromElectronsInPackedGen[1].e._eta<<" phi:"<<darkPhotonFromElectronsInPackedGen[1].e._phi<<endl;

    pkdGenEleCutApplied_N = count_if(cbegin(*pkdGenHdl_), cend(*pkdGenHdl_),
            [](const auto& p){ return abs(p.pdgId()) == 11 && p.et() > 4 && p.pt() > 2 && abs(p.eta())<1.479; });
    pfEp_N = count_if(cbegin(*pfHdl_), cend(*pfHdl_),
            [](const auto& p){ return abs(p.pdgId()) == 11; });
    pfElectron_N = count_if(cbegin(*pfHdl_), cend(*pfHdl_),
            [](const auto& p){ return abs(p.pdgId()) == 11 && p.charge()<0; });
    pfPositron_N = count_if(cbegin(*pfHdl_), cend(*pfHdl_),
            [](const auto& p){ return abs(p.pdgId()) == 11 && p.charge()>0; });
    pfGamma_N = count_if(cbegin(*pfHdl_), cend(*pfHdl_),
            [](const auto& p){ return p.pdgId() == 22; });
    patJet_N = patJetHdl_->size();
    darkPhotonWithAtLeastOneDaughtermatched_N = 0;
    matchedDarkPhoton_N = 0;

    if (darkPhotonFromElectronsInPackedGen.size() == 2) {

        // Collect edm ptrs for electrons/positrons which form darkphotons.
        vector<Ptr<pat::PackedGenParticle> > electronFromGenZp{}, positronFromGenZp{};
        for (const auto& q : darkPhotonFromElectronsInPackedGen) {
            electronFromGenZp.emplace_back(pkdGenHdl_, q.e.indexInCollection());
            positronFromGenZp.emplace_back(pkdGenHdl_, q.p.indexInCollection());
        }
        
        // Collect edm ptrs for electrons/positrons in PackedCandidate(reco).
        vector<Ptr<pat::PackedCandidate> >   electronInPfCands{}, positronInPfCands{};
        for (const auto& q : pfPtr_) {
            if (abs(q->pdgId()) != 11) continue;
            if (q->charge()<0) electronInPfCands.push_back(q);
            if (q->charge()>0) positronInPfCands.push_back(q);
        }

        ///DEBUG cout<<"Event"<<eventNum_<<"\t";
        if (electronInPfCands.size()>0 && electronFromGenZp.size()>0) {
            sidm::pairvec<pat::PackedCandidate, pat::PackedGenParticle> RecoGenElectron(electronInPfCands, electronFromGenZp);
            vector<pair<Ptr<pat::PackedCandidate>, Ptr<pat::PackedGenParticle> > > RecoGenElectronPair = RecoGenElectron.get();

            // Fill in the most matched (min dR) electron(s).
            sort(begin(RecoGenElectronPair), end(RecoGenElectronPair),
                    [](auto& lhs, auto& rhs) { return sidm::dR(lhs.first, lhs.second) < sidm::dR(rhs.first, rhs.second); });
            
            RecoGenElectron = sidm::pairvec<pat::PackedCandidate, pat::PackedGenParticle>(RecoGenElectronPair, true);
            RecoGenElectronPair.clear();
            RecoGenElectronPair = RecoGenElectron.get_zip();

            ///DEBUG cout<<"Ele  ";
            ///DEBUG cout<<electronInPfCands.size()<<" "<<electronFromGenZp.size()<<" ";
            ///DEBUG cout<<RecoGenElectronPair.size()<<" | ";
            for (const auto& q : RecoGenElectronPair) {
                electronFromPkdPF_ = sidm::Ep(q.first);
                mindRRecoWithGen = sidm::dR(q.first, q.second);
                ///DEBUG cout<<mindRRecoWithGen<<" ";
                electronsFromPackedPFTree_->Fill();

                for (auto& zp : darkPhotonFromElectronsInPackedGen) {
                    if (zp.e.indexInCollection() == q.second.key()) {
                        zp.e.matched = true;
                        zp.e.setIndexInCollectionMatched(q.first.key());
                    }
                }

            }
        }

        if (positronInPfCands.size()>0 && positronFromGenZp.size()>0) {
            sidm::pairvec<pat::PackedCandidate, pat::PackedGenParticle> RecoGenPositron(positronInPfCands, positronFromGenZp);
            vector<pair<Ptr<pat::PackedCandidate>, Ptr<pat::PackedGenParticle> > > RecoGenPositronPair = RecoGenPositron.get();

            // Fill in the most matched (min dR) positron(s).
            sort(begin(RecoGenPositronPair), end(RecoGenPositronPair),
                    [](const auto& lhs, const auto& rhs) { return sidm::dR(lhs.first, lhs.second) < sidm::dR(rhs.first, rhs.second); });

            RecoGenPositron = sidm::pairvec<pat::PackedCandidate, pat::PackedGenParticle>(RecoGenPositronPair, true);
            RecoGenPositronPair.clear();
            RecoGenPositronPair = RecoGenPositron.get_zip();

            ///DEBUG cout<<"Pos  ";
            ///DEBUG cout<<positronInPfCands.size()<<" "<<positronFromGenZp.size()<<" ";
            ///DEBUG cout<<RecoGenPositronPair.size()<<" | ";
            for (const auto& q : RecoGenPositronPair) {
                electronFromPkdPF_ = sidm::Ep(q.first);
                mindRRecoWithGen = sidm::dR(q.first, q.second);
                ///DEBUG cout<<mindRRecoWithGen<<" ";
                electronsFromPackedPFTree_->Fill();
                
                for (auto& zp : darkPhotonFromElectronsInPackedGen) {
                    if (zp.p.indexInCollection() == q.second.key()) {
                        zp.p.matched = true;
                        zp.p.setIndexInCollectionMatched(q.first.key());
                    }
                }

            }
        }
        ///DEBUG cout<<endl;

        matchedDarkPhoton_N = count_if(cbegin(darkPhotonFromElectronsInPackedGen), cend(darkPhotonFromElectronsInPackedGen),
                [](const auto& zp){return zp.e.matched && zp.p.matched;});
        darkPhotonWithAtLeastOneDaughtermatched_N = count_if(cbegin(darkPhotonFromElectronsInPackedGen), cend(darkPhotonFromElectronsInPackedGen),
                [](const auto& zp){return zp.e.matched || zp.p.matched;});
        matchedDarkPhotonWithJetIncluded_N = matchedDarkPhoton_N;


        for (const auto& zp : darkPhotonFromElectronsInPackedGen) {
            if (zp.e.matched && zp.p.matched) continue;

            /// Fill in electrons who missed reco.
            if (zp.e.matched == false) {
                genElectronNoReco_ = sidm::Ep(zp.e);
                genElectronNoRecoTree_->Fill();
            }
            if (zp.p.matched == false) {
                genElectronNoReco_ = sidm::Ep(zp.p);
                genElectronNoRecoTree_->Fill();
            }

            /// Count number of dark photons who has jet in the cone also.
            mindRZpWithJet = 9999.;
            int matched_jet_id = -1;
            for (const auto& j : patJetPtr_) {
                if ( sidm::dR(&zp, j) < mindRZpWithJet ) {
                    mindRZpWithJet = sidm::dR(&zp, j);
                    matched_jet_id = j.key();
                }
            }
            if ( mindRZpWithJet > 3. ) continue;

            Ptr<pat::Jet> matchedJetPtr(patJetHdl_, matched_jet_id);
            ++matchedDarkPhotonWithJetIncluded_N;
            suspiciousPatJet_ = sidm::Jet(matchedJetPtr);
            suspiciousPatJetTree_->Fill();
            cout<<"Event"<<eventNum_<<" min dR(zp, jet) = "<<mindRZpWithJet<<endl;
            if (!zp.e.matched && !zp.p.matched) {
                suspiciousPatJetSingle_ = sidm::Jet(matchedJetPtr);
                suspiciousPatJetSingleTree_->Fill();
            } else {
                suspiciousPatJetCoex_ = sidm::Jet(matchedJetPtr);
                //float dRJetWithMatchedEp = 9999.;
                if (zp.e.matched) {
                    dRJetWithMatchedEp = sidm::dR(matchedJetPtr, &zp.e);
                    for (auto& id : suspiciousPatJetCoex_._daughter_id) {
                        if (id == zp.e.indexInCollectionMatched()) {
                            ++epOfZpInJet_N;
                            break;
                        }
                    }
                } else {
                    dRJetWithMatchedEp = sidm::dR(matchedJetPtr, &zp.p);
                    for (auto& id : suspiciousPatJetCoex_._daughter_id) {
                        if (id == zp.p.indexInCollectionMatched()) {
                            ++epOfZpInJet_N;
                            break;
                        }
                    }
                }

                suspiciousPatJetCoexTree_->Fill();

            }

            //matchedDarkPhotonWithJetIncluded_N +=
            //    count_if(cbegin(patJetPtr_), cend(patJetPtr_), [zp](const auto& j){ return sidm::dR(&zp, j)<=0.4; });
        }
        //cout<<"Event"<<eventNum_<<": matched darkphoton "<<matchedDarkPhoton_N<<" "<<darkPhotonWithAtLeastOneDaughtermatched_N
        //    <<" "<<matchedDarkPhotonWithJetIncluded_N<<endl;


        vector<Ptr<pat::Electron> > patElectronPtr_ = patElectronHdl_->ptrs();
        for (const auto& p : patElectronPtr_) {
            electronFromPat_ = sidm::Ep(p);

            vector<float> deltaRtmp{};
            if (p->charge()<0) {
                for (const auto& zp : darkPhotonFromElectronsInPackedGen) {
                    deltaRtmp.push_back(sidm::dR(p, &zp.e));
                }
            } else if (p->charge()>0) {
                for (const auto& zp : darkPhotonFromElectronsInPackedGen) {
                    deltaRtmp.push_back(sidm::dR(p, &zp.p));
                }
            }
            sort(deltaRtmp.begin(), deltaRtmp.end());
            mindRPatWithGen = deltaRtmp[0];
            if (mindRPatWithGen>15) continue;

            electronsFromPATTree_->Fill();
        }
    } else {
        cout<<"Event"<<eventNum_<<": NOT 2 darkphotons are reconstructed from GEN, instead-  "
            <<darkPhotonFromElectronsInPackedGen.size()<<".\n";
    }

    eventTree_->Fill();
    ++eventNum_;
}


// ------------ method called once each job just before starting event loop  ------------
void 
sidm::electronFinder::beginJob()
{
    /// Statistics of multiplicities of collections per event
    eventTree_ = fs_->make<TTree>("eventTree", "information per event");
    
    eventTree_->Branch("eventId", &eventNum_, "eventId/I");
    eventTree_->Branch("numOfEleInPackedGen", &pkdGenElectron_N, "numOfEleInPackedGen/I");
    eventTree_->Branch("numOfEleInPackedGenWithCut", &pkdGenEleCutApplied_N, "numOfEleInPackedGenWithCut/I");
    eventTree_->Branch("numOfEpInPfCands",  &pfEp_N, "numOfEpInPfCands/I");
    eventTree_->Branch("numOfEleInPfCands", &pfElectron_N, "numOfEleInPfCands/I");
    eventTree_->Branch("numOfPosInPfCands", &pfPositron_N, "numOfPosInPfCands/I");
    eventTree_->Branch("numOfMatchedDarkPhotons", &matchedDarkPhoton_N, "numOfMatchedDarkPhotons/I");
    eventTree_->Branch("numOfDarkPhotonsWithAtLeastOneDaughterMatched", &darkPhotonWithAtLeastOneDaughtermatched_N, "numOfDarkPhotonsWithAtLeastOneDaughterMatched/I");
    eventTree_->Branch("numOfJets", &patJet_N, "numOfJets/I");
    eventTree_->Branch("numOfMatchedDarkPhotonWithJetIncluded", &matchedDarkPhotonWithJetIncluded_N, "numOfMatchedDarkPhotonWithJetIncluded/I");
    //eventTree_->Branch("numOfPhoInPfCands", &pfGamma_N, "numOfPhoInPfCands/I");

    darkPhotonFromGenElectronsTree_ = fs_->make<TTree>("darkPhotonFromGenElectrons",
            "darkphotons reco from electrons in packedGenParticles collection");
    darkPhotonFromGenElectronsTree_->Branch("eventId", &darkPhotonFromGenEle_._eventId, "eventId/I");
    darkPhotonFromGenElectronsTree_->Branch("pt",      &darkPhotonFromGenEle_._pt, "pt/F");
    darkPhotonFromGenElectronsTree_->Branch("eta",     &darkPhotonFromGenEle_._eta, "eta/F");
    darkPhotonFromGenElectronsTree_->Branch("mass",    &darkPhotonFromGenEle_._mass, "mass/F");
    darkPhotonFromGenElectronsTree_->Branch("invm",    &darkPhotonFromGenEle_._invM, "invm/F");
    darkPhotonFromGenElectronsTree_->Branch("pt_e",    &darkPhotonFromGenEle_.e._pt, "pt_e/F");
    darkPhotonFromGenElectronsTree_->Branch("pt_p",    &darkPhotonFromGenEle_.p._pt, "pt_p/F");
    darkPhotonFromGenElectronsTree_->Branch("et_e",    &darkPhotonFromGenEle_.e._et, "et_e/F");
    darkPhotonFromGenElectronsTree_->Branch("et_p",    &darkPhotonFromGenEle_.p._et, "et_p/F");
    darkPhotonFromGenElectronsTree_->Branch("eta_e",   &darkPhotonFromGenEle_.e._eta, "eta_e/F");
    darkPhotonFromGenElectronsTree_->Branch("eta_p",   &darkPhotonFromGenEle_.p._eta, "eta_p/F");
    darkPhotonFromGenElectronsTree_->Branch("dR_ep",   &darkPhotonFromGenEle_._dR, "dR_ep/F");
    darkPhotonFromGenElectronsTree_->Branch("dEta_ep", &darkPhotonFromGenEle_._dEta, "dEta_ep/F");
    darkPhotonFromGenElectronsTree_->Branch("dPhi_ep", &darkPhotonFromGenEle_._dPhi, "dPhi_ep/F");
    darkPhotonFromGenElectronsTree_->Branch("dv_x",    &darkPhotonFromGenEle_._dv_x, "dv_x/D");
    darkPhotonFromGenElectronsTree_->Branch("dv_y",    &darkPhotonFromGenEle_._dv_y, "dv_y/D");
    darkPhotonFromGenElectronsTree_->Branch("dv_z",    &darkPhotonFromGenEle_._dv_z, "dv_z/D");

    electronsFromPackedPFTree_ = fs_->make<TTree>("electronsFromPackedPF", "electron info from packedPF collection");
    electronsFromPackedPFTree_->Branch("eventId", &eventNum_, "eventId/I");
    electronsFromPackedPFTree_->Branch("pt", &electronFromPkdPF_._pt, "pt/F");
    electronsFromPackedPFTree_->Branch("eta", &electronFromPkdPF_._eta, "eta/F");
    electronsFromPackedPFTree_->Branch("et", &electronFromPkdPF_._et, "et/F");
    electronsFromPackedPFTree_->Branch("dR", &mindRRecoWithGen, "dR/F");

    electronsFromPATTree_ = fs_->make<TTree>("electronsFromPAT", "electron info from pat::Electron collection");
    electronsFromPATTree_->Branch("eventId", &eventNum_, "eventId/I");
    electronsFromPATTree_->Branch("pt", &electronFromPat_._pt, "pt/F");
    electronsFromPATTree_->Branch("eta", &electronFromPat_._eta, "eta/F");
    electronsFromPATTree_->Branch("et", &electronFromPat_._et, "et/F");
    electronsFromPATTree_->Branch("dR", &mindRPatWithGen, "dR/F");

    genElectronNoRecoTree_ = fs_->make<TTree>("genElectronNoReco", "electron info from packedGenParticle collection\
            who missed the reconstruction");
    genElectronNoRecoTree_->Branch("eventId", &eventNum_, "eventId/I");
    genElectronNoRecoTree_->Branch("pt",  &genElectronNoReco_._pt, "pt/F");
    genElectronNoRecoTree_->Branch("eta", &genElectronNoReco_._eta, "eta/F");
    genElectronNoRecoTree_->Branch("et",  &genElectronNoReco_._et, "et/F");
    
    suspiciousPatJetTree_ = fs_->make<TTree>("suspiciousPatJet", "jet info who is inside dark photon cone");
    suspiciousPatJetTree_->Branch("eventId", &eventNum_, "eventId/I");
    suspiciousPatJetTree_->Branch("pt",  &suspiciousPatJet_._pt, "pt/F");
    suspiciousPatJetTree_->Branch("eta", &suspiciousPatJet_._eta, "eta/F");
    suspiciousPatJetTree_->Branch("et",  &suspiciousPatJet_._et, "et/F");
    suspiciousPatJetTree_->Branch("chargedHOverE",  &suspiciousPatJet_._charged_H_over_E, "chargedHOverE/F");
    suspiciousPatJetTree_->Branch("chargedMultiplicity",  &suspiciousPatJet_._chargedMultiplicity, "chargedMultiplicity/I");
    suspiciousPatJetTree_->Branch("electronEnergyFraction",  &suspiciousPatJet_._electronEnergyFraction, "electronEnergyFraction/F");
    suspiciousPatJetTree_->Branch("electronMultiplicity",  &suspiciousPatJet_._electronMultiplicity, "electronMultiplicity/I");
    suspiciousPatJetTree_->Branch("numberOfDaughters",  &suspiciousPatJet_._num_of_daughters, "numberOfDaughers/I");
    suspiciousPatJetTree_->Branch("daughterIds",  "vector<unsigned int>", &suspiciousPatJet_._daughter_id);
    suspiciousPatJetTree_->Branch("dRZpWithJet",  &mindRZpWithJet, "dRZpWithJet/F");

    suspiciousPatJetSingleTree_ = fs_->make<TTree>("suspiciousPatJetSingle", "jet info who is inside dark photon cone with no electron or positron along");
    suspiciousPatJetSingleTree_->Branch("eventId", &eventNum_, "eventId/I");
    suspiciousPatJetSingleTree_->Branch("pt",  &suspiciousPatJetSingle_._pt, "pt/F");
    suspiciousPatJetSingleTree_->Branch("eta", &suspiciousPatJetSingle_._eta, "eta/F");
    suspiciousPatJetSingleTree_->Branch("et",  &suspiciousPatJetSingle_._et, "et/F");
    suspiciousPatJetSingleTree_->Branch("chargedHOverE",  &suspiciousPatJetSingle_._charged_H_over_E, "chargedHOverE/F");
    suspiciousPatJetSingleTree_->Branch("chargedMultiplicity",  &suspiciousPatJetSingle_._chargedMultiplicity, "chargedMultiplicity/I");
    suspiciousPatJetSingleTree_->Branch("electronEnergyFraction",  &suspiciousPatJetSingle_._electronEnergyFraction, "electronEnergyFraction/F");
    suspiciousPatJetSingleTree_->Branch("electronMultiplicity",  &suspiciousPatJetSingle_._electronMultiplicity, "electronMultiplicity/I");
    suspiciousPatJetSingleTree_->Branch("numberOfDaughters",  &suspiciousPatJetSingle_._num_of_daughters, "numberOfDaughers/I");
    suspiciousPatJetSingleTree_->Branch("daughterIds",  "vector<unsigned int>", &suspiciousPatJetSingle_._daughter_id);

    suspiciousPatJetCoexTree_ = fs_->make<TTree>("suspiciousPatJetCoex", "jet info who is inside dark photon cone with electron or positron along");
    suspiciousPatJetCoexTree_->Branch("eventId", &eventNum_, "eventId/I");
    suspiciousPatJetCoexTree_->Branch("pt",  &suspiciousPatJetCoex_._pt, "pt/F");
    suspiciousPatJetCoexTree_->Branch("eta", &suspiciousPatJetCoex_._eta, "eta/F");
    suspiciousPatJetCoexTree_->Branch("et",  &suspiciousPatJetCoex_._et, "et/F");
    suspiciousPatJetCoexTree_->Branch("chargedHOverE",  &suspiciousPatJetCoex_._charged_H_over_E, "chargedHOverE/F");
    suspiciousPatJetCoexTree_->Branch("chargedMultiplicity",  &suspiciousPatJetCoex_._chargedMultiplicity, "chargedMultiplicity/I");
    suspiciousPatJetCoexTree_->Branch("electronEnergyFraction",  &suspiciousPatJetCoex_._electronEnergyFraction, "electronEnergyFraction/F");
    suspiciousPatJetCoexTree_->Branch("electronMultiplicity",  &suspiciousPatJetCoex_._electronMultiplicity, "electronMultiplicity/I");
    suspiciousPatJetCoexTree_->Branch("numberOfDaughters",  &suspiciousPatJetCoex_._num_of_daughters, "numberOfDaughers/I");
    suspiciousPatJetCoexTree_->Branch("daughterIds",  "vector<unsigned int>", &suspiciousPatJetCoex_._daughter_id);
    suspiciousPatJetCoexTree_->Branch("JetDaughterIsMatchedEp", &epOfZpInJet_N, "JetDaughterIsMatchedEp/I");
    suspiciousPatJetCoexTree_->Branch("dRJetMatchdEp", &dRJetWithMatchedEp, "dRJetMatchdEp/F");

}

// ------------ method called once each job just after ending the event loop  ------------
void 
sidm::electronFinder::endJob() 
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
sidm::electronFinder::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(sidm::electronFinder);
