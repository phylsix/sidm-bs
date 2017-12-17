// -*- C++ -*-
//
// Package:    sidm-bs/MuonAnalysis
// Class:      MuonAnalysis
#include <algorithm>
#include <cmath>
#include <map>

#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/Common/interface/Ptr.h"
#include "DataFormats/Math/interface/LorentzVector.h"
#include "DataFormats/Math/interface/deltaR.h"

#include "sidm-bs/jetComponent/interface/physicsObject.h"
#include "sidm-bs/jetComponent/interface/utilities.h"

#include "sidm-bs/jetComponent/interface/MuonAnalysis.h"

sidm::MuonAnalysis::MuonAnalysis(const edm::ParameterSet& iConfig):
    genParticleTk_(consumes<edm::View<reco::GenParticle> >(iConfig.getUntrackedParameter<edm::InputTag>("GenParticleTag_", edm::InputTag("prunedGenParticles")))),
    pkdGenTk_(consumes<edm::View<pat::PackedGenParticle> >(iConfig.getUntrackedParameter<edm::InputTag>("PatPkdGenTag_", edm::InputTag("packedGenParticles")))),
    pfTk_(consumes<edm::View<pat::PackedCandidate> >(iConfig.getUntrackedParameter<edm::InputTag>("PatPfTag_", edm::InputTag("packedPFCandidates")))),
    patElectronTk_(consumes<edm::View<pat::Electron> >(iConfig.getUntrackedParameter<edm::InputTag>("PatElectronTag_", edm::InputTag("slimmedElectrons")))),
    patMuonTk_(consumes<edm::View<pat::Muon> >(iConfig.getUntrackedParameter<edm::InputTag>("PatMuonTag_", edm::InputTag("slimmedMuons")))),
    patJetTk_(consumes<edm::View<pat::Jet> >(iConfig.getUntrackedParameter<edm::InputTag>("PatJetTag_", edm::InputTag("slimmedJets"))))
{
    usesResource("TFileService");
    eventNum_ = 0;
}


sidm::MuonAnalysis::~MuonAnalysis()
{

    std::cout << "\n\nNumber of events: " << eventNum_ << "\n\n";
    eventNum_ = 0;

}

// ------------ method called for each event  ------------
void
sidm::MuonAnalysis::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
    using namespace edm;
    using namespace std;

    Handle<View<reco::GenParticle> > genParticleHdl_;
    iEvent.getByToken(genParticleTk_, genParticleHdl_);

    Handle<View<pat::PackedGenParticle> > pkdGenHdl_;
    iEvent.getByToken(pkdGenTk_, pkdGenHdl_);

    Handle<View<pat::Electron> > patElectronHdl_;
    iEvent.getByToken(patElectronTk_, patElectronHdl_);

    Handle<View<pat::Muon> > patMuonHdl_;
    iEvent.getByToken(patMuonTk_, patMuonHdl_);

    Handle<View<pat::Jet> > patJetHdl_;
    iEvent.getByToken(patJetTk_, patJetHdl_);

    Handle<View<pat::PackedCandidate> > pfHdl_;
    iEvent.getByToken(pfTk_, pfHdl_);

    vector<Ptr<pat::Muon> >             patMuonPtr_ = patMuonHdl_->ptrs();
    vector<Ptr<pat::Jet> >               patJetPtr_ = patJetHdl_->ptrs();
    vector<Ptr<reco::GenParticle> >      genPtr_    = genParticleHdl_->ptrs();
    vector<Ptr<pat::PackedGenParticle> > pkdGenPtr_ = pkdGenHdl_->ptrs();

    //>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> PAT muon >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

    muonMultiplicity_->Fill(patMuonHdl_->size());
    for (const auto& _mu : *patMuonHdl_) {
        muonEnergy_->Fill(_mu.energy());
        muonPt_->Fill(_mu.pt());
        muonEt_->Fill(_mu.et());
        muonEta_->Fill(_mu.eta());
        muonPhi_->Fill(_mu.phi());
        muonVertex_->Fill(abs(_mu.vz()), sqrt(_mu.vx()*_mu.vx() + _mu.vy()*_mu.vy()));
    }

    //<<<<<<<<<<<<<<<<<<<<<<<<<<<<< Gen Particle <<<<<<<<<<<<<<<<<<<<<<<<<<<<< 

    //>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> Gen Particle >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

    /**
     * _genDarkphotos_
     * KEY:   The key of darkPhoton in reco::GenParticle collection
     * VALUE: Pair of keys of corresponding muon daughters in pat::PackedGenParticles collection,
     *        first key is muon+, second is muon-.
     * 
     * Keys can be used to uniquely construct edm::Ptr together with Handle so we
     * can keep track.
     */
    map<unsigned int, pair<unsigned int, unsigned int> > genDarkphotons;
    genDarkphotons.clear();

    for (const auto& _dp : genPtr_) {
        if (_dp->pdgId() != 32) continue;
        vector<Ptr<pat::PackedGenParticle> > legitMuonDaughters{};
        
        for (const auto& _gMu : pkdGenPtr_) {
            if (_gMu->status()!=1 || abs(_gMu->pdgId())!=13) continue;
            const reco::Candidate* motherInPrunedCollection(_gMu->mother(0));
            if (motherInPrunedCollection!=nullptr && sidm::is_ancestor(_dp.get(),motherInPrunedCollection)) {
                legitMuonDaughters.push_back(_gMu);
            }
        }

        if (legitMuonDaughters.size()!=2) {
            cout<<"Event"<<setw(5)<<left<<eventNum_<<__LINE__<<"Dark photon has "<<legitMuonDaughters.size()<<" muon daughters instead of 2.\n";
            continue;
        }

        if (legitMuonDaughters[0]->charge()*legitMuonDaughters[1]->charge() != -1) {
            cout<<"Event"<<setw(5)<<left<<eventNum_<<__LINE__<<"Dark photon has same sign muon daughters!\n";
            continue;
        }

        if (legitMuonDaughters[0]->charge() > 0) {
           genDarkphotons[_dp.key()] = std::make_pair(legitMuonDaughters[0].key(), legitMuonDaughters[1].key());
        } else {
           genDarkphotons[_dp.key()] = std::make_pair(legitMuonDaughters[1].key(), legitMuonDaughters[0].key());
        }

        for (const auto& _gMu : legitMuonDaughters) {
            genMuonEnergy_->Fill(_gMu->energy());
            genMuonPt_    ->Fill(_gMu->pt());
            genMuonEt_    ->Fill(_gMu->et());
            genMuonEta_   ->Fill(_gMu->eta());
            genMuonPhi_   ->Fill(_gMu->phi());
            genMuonVertex_->Fill(abs(_dp->daughter(0)->vz()),
                                sqrt(_dp->daughter(0)->vx()*_dp->daughter(0)->vx() +
                                      _dp->daughter(0)->vy()*_dp->daughter(0)->vy()));
            genMuonVertex_->Fill(abs(_dp->daughter(0)->vz()),
                                sqrt(_dp->daughter(0)->vx()*_dp->daughter(0)->vx() +
                                      _dp->daughter(0)->vy()*_dp->daughter(0)->vy()));

        }

    }
    genMuonMultiplicity_->Fill(2*genDarkphotons.size());

    if (genDarkphotons.size() != 2) {
        cout<<"Event"<<setw(5)<<left<<eventNum_<<__LINE__<<genDarkphotons.size()<<" gen dark photons found instead of 2.\n";
        return;
    }

    //<<<<<<<<<<<<<<<<<<<<<<<<<<<<< Gen Particle <<<<<<<<<<<<<<<<<<<<<<<<<<<<< 



    //>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> Gen PAT Matching >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    map<unsigned int, pair<bool, bool> > genDarkphotonsMatchFlags;
    vector<Ptr<pat::PackedGenParticle> > genMuPlus{}, genMuMinus{};
    for (const auto& gd : genDarkphotons) {
        genDarkphotonsMatchFlags[gd.first] = make_pair(false, false);
        
        genMuPlus.emplace_back(pkdGenHdl_, gd.second.first);
        genMuMinus.emplace_back(pkdGenHdl_, gd.second.second);
    }

    vector<Ptr<pat::Muon> > patMuPlus{}, patMuMinus{};
    for (const auto& mu : patMuonPtr_) {
        if (mu->charge()<0) patMuPlus.push_back(mu);
        else patMuMinus.push_back(mu);
    }


    /**
     * _MuonPatGenMap_
     * KEY:   key in pat::Muon collection
     * VALUE: key in pat::PackedGenParticle collection.
     */
    map<unsigned int, unsigned int> MuonPatGenMap;
    
    /*Mu plus*/
    if (genMuPlus.size()>0 && patMuPlus.size()>0) {
        sidm::pairvec<pat::Muon, pat::PackedGenParticle> PatGenMuPlus(patMuPlus, genMuPlus);
        vector<pair<Ptr<pat::Muon>, Ptr<pat::PackedGenParticle> > > PatGenMuPlusPairs = PatGenMuPlus.get();
        sort(begin(PatGenMuPlusPairs), end(PatGenMuPlusPairs), [](auto& lhs, auto& rhs)
                {return sidm::dR(lhs.first,lhs.second) < sidm::dR(rhs.first,rhs.second);});

        PatGenMuPlus = sidm::pairvec<pat::Muon, pat::PackedGenParticle>(PatGenMuPlusPairs,true);
        PatGenMuPlusPairs.clear();
        PatGenMuPlusPairs = PatGenMuPlus.get_zip();

        for (const auto& pg : PatGenMuPlusPairs) {
            MuonPatGenMap[pg.first.key()] = pg.second.key();
            /*Fill the matched pat muon properties+dR*/
            matchedMuonEnergy_->Fill(pg.first->energy());
            matchedMuonPt_    ->Fill(pg.first->pt());
            matchedMuonEt_    ->Fill(pg.first->et());
            matchedMuonEta_   ->Fill(pg.first->eta());
            matchedMuonPhi_   ->Fill(pg.first->phi());
            matchedMuonVertex_->Fill(abs(pg.first->vz()),
                                    sqrt(pg.first->vx()*pg.first->vx() + pg.first->vy()*pg.first->vy()));
            matchedMuonDR_    ->Fill(sidm::dR(pg.first, pg.second));
            /*---------------------------------------*/
            for (const auto& gd : genDarkphotons) {
                if (gd.second.first == pg.second.key()) genDarkphotonsMatchFlags[gd.first].first = true;
            }
        }
    }

    /*Mu minus*/
    if (genMuMinus.size()>0 && patMuMinus.size()>0) {
        //assert(0);
        sidm::pairvec<pat::Muon, pat::PackedGenParticle> PatGenMuMinus(patMuMinus, genMuMinus);
        vector<pair<Ptr<pat::Muon>, Ptr<pat::PackedGenParticle> > > PatGenMuMinusPairs = PatGenMuMinus.get();
        sort(begin(PatGenMuMinusPairs), end(PatGenMuMinusPairs), [](auto& lhs, auto& rhs)
                {return sidm::dR(lhs.first,lhs.second) < sidm::dR(rhs.first,rhs.second);});

        PatGenMuMinus = sidm::pairvec<pat::Muon, pat::PackedGenParticle>(PatGenMuMinusPairs,true);
        PatGenMuMinusPairs.clear();
        PatGenMuMinusPairs = PatGenMuMinus.get_zip();

        for (const auto& pg : PatGenMuMinusPairs) {
            MuonPatGenMap[pg.first.key()] = pg.second.key();
            /*Fill the matched pat muon properties+dR*/
            matchedMuonEnergy_->Fill(pg.first->energy());
            matchedMuonPt_    ->Fill(pg.first->pt());
            matchedMuonEt_    ->Fill(pg.first->et());
            matchedMuonEta_   ->Fill(pg.first->eta());
            matchedMuonPhi_   ->Fill(pg.first->phi());
            matchedMuonVertex_->Fill(abs(pg.first->vz()),
                                    sqrt(pg.first->vx()*pg.first->vx() + pg.first->vy()*pg.first->vy()));
            matchedMuonDR_    ->Fill(sidm::dR(pg.first, pg.second));
            /*---------------------------------------*/
            for (const auto& gd : genDarkphotons) {
                if (gd.second.second == pg.second.key()) genDarkphotonsMatchFlags[gd.first].second = true;
            }
        }
    }

    int totalMatchedMuons(0);
    for (const auto& gd : genDarkphotonsMatchFlags) {
        if (gd.second.first) ++totalMatchedMuons;
        if (gd.second.second) ++totalMatchedMuons;
    }
    assert(totalMatchedMuons<=4);
    matchedMuonMultiplicity_->Fill(totalMatchedMuons);

    //<<<<<<<<<<<<<<<<<<<<<<<<<<<<< Gen PAT Matching <<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
    


    //>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> Unmatched Gen Muons >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    
    for (const auto& gd : genDarkphotonsMatchFlags) {
        /*Unmatched mu plus*/
        if (!gd.second.first) {
            Ptr<pat::PackedGenParticle> unMuPlus(pkdGenHdl_,genDarkphotons[gd.first].first);
            /*Fill unmatched gen muon properties*/
            unmatchedMuonEnergy_->Fill(unMuPlus->energy());
            unmatchedMuonPt_    ->Fill(unMuPlus->pt());
            unmatchedMuonEt_    ->Fill(unMuPlus->et());
            unmatchedMuonEta_   ->Fill(unMuPlus->eta());
            unmatchedMuonPhi_   ->Fill(unMuPlus->phi());
            unmatchedMuonVertex_->Fill(abs(unMuPlus->vz()),
                                      sqrt(unMuPlus->vx()*unMuPlus->vx() + unMuPlus->vy()*unMuPlus->vy()));
        }

        /*Unmatched mu minus*/
        if (!gd.second.second) {
            Ptr<pat::PackedGenParticle> unMuMinus(pkdGenHdl_, genDarkphotons[gd.first].second);
            /*Fill unmatched gen muon properties*/
            unmatchedMuonEnergy_->Fill(unMuMinus->energy());
            unmatchedMuonPt_    ->Fill(unMuMinus->pt());
            unmatchedMuonEt_    ->Fill(unMuMinus->et());
            unmatchedMuonEta_   ->Fill(unMuMinus->eta());
            unmatchedMuonPhi_   ->Fill(unMuMinus->phi());
            unmatchedMuonVertex_->Fill(abs(unMuMinus->vz()),
                                  sqrt(unMuMinus->vx()*unMuMinus->vx() + unMuMinus->vy()*unMuMinus->vy()));
        }
    }
    unmatchedMuonMultiplicity_->Fill(4-totalMatchedMuons);
    //<<<<<<<<<<<<<<<<<<<<<<<<<<<<< Unmatched Gen Muons <<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
    //
    //
    //
    //>>>>>>>>>>>>>>>>>>>>>>>>>>>>> Matched PAT Muon Pair >>>>>>>>>>>>>>>>>>>>>>>>>>>>>

    int totalMatchedDarkphotons(0);
    for (const auto& gd : genDarkphotonsMatchFlags) {
        if (gd.second.first && gd.second.second) {
            ++totalMatchedDarkphotons;
            /*Fill invM + dR*/
            Ptr<pat::PackedGenParticle> muPlus(pkdGenHdl_, genDarkphotons[gd.first].first);
            Ptr<pat::PackedGenParticle> muMinus(pkdGenHdl_, genDarkphotons[gd.first].second);

            matchedDarkphotonInvm_MuMu_->Fill((muPlus->p4()+muMinus->p4()).M());
            matchedDarkphotonDR_MuMu_  ->Fill(sidm::dR(muPlus, muMinus));
        }
    }
    assert(totalMatchedDarkphotons<=2);
    matchedDarkphotonMultiplicity_MuMu_->Fill(totalMatchedDarkphotons);
    //<<<<<<<<<<<<<<<<<<<<<<<<<<<<< Matched PAT Muon Pair <<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
    //
    //
    //
    //>>>>>>>>>>>>>>>>>>>>>>>>>>>>> Matched PAT Muon Pair >>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    
    map<unsigned int, unsigned int> genDarkphotonJetInConeMap;
    int aloneJetMatched(0);
    int matchedDarkphoton_MuJet(0);

    for (auto& gd : genDarkphotonsMatchFlags) {
        if (gd.second.first && gd.second.second) { continue; }
        /* Check if there is a jet within a cone of dR<3, and pick the closest one */
        float darkphotonJetDeltaR_Min(9999.);
        unsigned int darkphotonClosestJetId(-1);
        Ptr<reco::GenParticle> genDarkphotonPtr(genParticleHdl_, gd.first);
        for (const auto& j : patJetPtr_) {
            float darkphotonJetDeltaR = sidm::dR(genDarkphotonPtr, j);
            if (darkphotonJetDeltaR < darkphotonJetDeltaR_Min) {
                darkphotonJetDeltaR_Min = darkphotonJetDeltaR;
                darkphotonClosestJetId = j.key();
            }
        }
        if (darkphotonJetDeltaR_Min > 3.) { continue; } // Too far away
        genDarkphotonJetInConeMap[gd.first] = darkphotonClosestJetId;
        Ptr<pat::Jet> closestJetPtr(patJetHdl_, darkphotonClosestJetId);
        /* Fill darkphotonJetDeltaR_Min + jet properties*/
        jetInConeEnergy_              ->Fill(closestJetPtr->energy());
        jetInConePt_                  ->Fill(closestJetPtr->pt()); 
        jetInConeEt_                  ->Fill(closestJetPtr->et());
        jetInConeEta_                 ->Fill(closestJetPtr->eta());
        jetInConePhi_                 ->Fill(closestJetPtr->phi());
        jetInConeChargedHoverE_       ->Fill( (closestJetPtr->isPFJet()||closestJetPtr->isJPTJet()) ? closestJetPtr->chargedHadronEnergy()/closestJetPtr->chargedEmEnergy() : -1 );
        jetInConeChargedMultiplicity_ ->Fill( (closestJetPtr->isPFJet()||closestJetPtr->isJPTJet()) ? closestJetPtr->chargedMultiplicity() : -1 );
        jetInConeMuonEnergyFraction_  ->Fill( closestJetPtr->isPFJet() ? closestJetPtr->muonEnergyFraction() : -1);
        jetInConeMuonMultiplicity_    ->Fill( closestJetPtr->isPFJet() ? closestJetPtr->muonMultiplicity()   : -1);
        jetInConeDR_                  ->Fill(darkphotonJetDeltaR_Min);
        /*---------------------------------------------*/
        if (gd.second.first==false && gd.second.second==false) { ++aloneJetMatched; }
        else {
            if (gd.second.first==false) {
                for (unsigned int i=0; i<closestJetPtr->numberOfDaughters(); ++i) {
                    if (closestJetPtr->daughter(i) != nullptr &&
                        closestJetPtr->daughter(i)->charge() > 0 &&
                        sidm::is_ancestor(genDarkphotonPtr.get(), closestJetPtr->daughter(i)))
                    {
                        ++matchedDarkphoton_MuJet;
                        break;
                    }
                }
            }
            if (gd.second.second==false) {
                for (unsigned int i=0; i<closestJetPtr->numberOfDaughters(); ++i) {
                    if (closestJetPtr->daughter(i) != nullptr &&
                        closestJetPtr->daughter(i)->charge() < 0 &&
                        sidm::is_ancestor(genDarkphotonPtr.get(), closestJetPtr->daughter(i)))
                    {
                        ++matchedDarkphoton_MuJet;
                        break;
                    }
                }
            }
        }
    }
    jetInConeMultiplicity_            ->Fill(genDarkphotonJetInConeMap.size());
    if (aloneJetMatched>0)         { jetInConeAloneMultiplicity_         ->Fill(aloneJetMatched); }
    if (matchedDarkphoton_MuJet>0) { matchedDarkphotonMultiplicity_MuJet_->Fill(matchedDarkphoton_MuJet); }
    ++eventNum_;
}


// ------------ method called once each job just before starting event loop  ------------
void 
sidm::MuonAnalysis::beginJob()
{
    /// Statistics of multiplicities of collections per event
    eventTree_ = fs_->make<TTree>("eventTree", "information per event");
    
    TFileDirectory muonDir = fs_->mkdir("PAT_Muon");
    muonMultiplicity_ = muonDir.make<TH1F>("MuonMultiplicity", "Muon Multiplicity", 10,0,10);
    muonEnergy_       = muonDir.make<TH1F>("MuonEnergy", "Muon Energy", 100,0,100);
    muonPt_           = muonDir.make<TH1F>("MuonPt", "Muon Pt", 100,0,100);
    muonEt_           = muonDir.make<TH1F>("MuonEt", "Muon Et", 100,0,100);
    muonEta_          = muonDir.make<TH1F>("MuonEta", "Muon Eta",60,-3,3);
    muonPhi_          = muonDir.make<TH1F>("MuonPhi", "Muon Phi",50,-3.1416,3.1416);
    muonVertex_       = muonDir.make<TH2F>("MuonVertex", "Muon vertex position",100,0,50,50,0,1);

    TFileDirectory genMuonDir = fs_->mkdir("GEN_Muon");
    genMuonMultiplicity_ = genMuonDir.make<TH1F>("genMuonMultiplicity", "Muon Multiplicity", 10,0,10);
    genMuonEnergy_       = genMuonDir.make<TH1F>("genMuonEnergy", "Muon Energy", 100,0,100);
    genMuonPt_           = genMuonDir.make<TH1F>("genMuonPt", "Muon Pt", 100,0,100);
    genMuonEt_           = genMuonDir.make<TH1F>("genMuonEt", "Muon Et", 100,0,100);
    genMuonEta_          = genMuonDir.make<TH1F>("genMuonEta", "Muon Eta",60,-3,3);
    genMuonPhi_          = genMuonDir.make<TH1F>("genMuonPhi", "Muon Phi",50,-3.1416,3.1416);
    genMuonVertex_       = genMuonDir.make<TH2F>("genMuonVertex", "Muon vertex position",100,0,50,50,0,1);

    TFileDirectory matchedMuonDir = fs_->mkdir("Matched_PAT_Muon");
    matchedMuonMultiplicity_ = matchedMuonDir.make<TH1F>("matchedMuonMultiplicity", "Muon Multiplicity", 10,0,10);
    matchedMuonEnergy_       = matchedMuonDir.make<TH1F>("matchedMuonEnergy", "Muon Energy", 100,0,100);
    matchedMuonPt_           = matchedMuonDir.make<TH1F>("matchedMuonPt", "Muon Pt", 100,0,100);
    matchedMuonEt_           = matchedMuonDir.make<TH1F>("matchedMuonEt", "Muon Et", 100,0,100);
    matchedMuonEta_          = matchedMuonDir.make<TH1F>("matchedMuonEta", "Muon Eta",60,-3,3);
    matchedMuonPhi_          = matchedMuonDir.make<TH1F>("matchedMuonPhi", "Muon Phi",50,-3.1416,3.1416);
    matchedMuonVertex_       = matchedMuonDir.make<TH2F>("matchedMuonVertex", "Muon vertex position",100,0,50,50,0,1);
    matchedMuonDR_           = matchedMuonDir.make<TH1F>("matchedMuonDR", "dR between pat muon and matched gen muon",100,0,0.5);

    TFileDirectory unmatchedMuonDir = fs_->mkdir("Unmatched_GEN_Muon");
    unmatchedMuonMultiplicity_ = unmatchedMuonDir.make<TH1F>("unmatchedMuonMultiplicity", "Muon Multiplicity", 10,0,10);
    unmatchedMuonEnergy_       = unmatchedMuonDir.make<TH1F>("unmatchedMuonEnergy", "Muon Energy", 100,0,100);
    unmatchedMuonPt_           = unmatchedMuonDir.make<TH1F>("unmatchedMuonPt", "Muon Pt", 100,0,100);
    unmatchedMuonEt_           = unmatchedMuonDir.make<TH1F>("unmatchedMuonEt", "Muon Et", 100,0,100);
    unmatchedMuonEta_          = unmatchedMuonDir.make<TH1F>("unmatchedMuonEta", "Muon Eta",60,-3,3);
    unmatchedMuonPhi_          = unmatchedMuonDir.make<TH1F>("unmatchedMuonPhi", "Muon Phi",50,-3.1416,3.1416);
    unmatchedMuonVertex_       = unmatchedMuonDir.make<TH2F>("unmatchedMuonVertex", "Muon vertex position",100,0,50,50,0,1);

    TFileDirectory matchedDarkphotonDir = fs_->mkdir("matched_darkphoton");
    matchedDarkphotonMultiplicity_MuMu_ = matchedDarkphotonDir.make<TH1F>("matchedDarkphotonMultiplicityMuMu", "Darkphoton(MuMu) Multiplicity",10,0,10);
    matchedDarkphotonInvm_MuMu_         = matchedDarkphotonDir.make<TH1F>("matchedDarkphotonInvmMuMu", "Darkphoton(MuMu) Invariant Mass",100,0,5);
    matchedDarkphotonDR_MuMu_           = matchedDarkphotonDir.make<TH1F>("matchedDarkphotonDRMuMu", "Darkphoton(MuMu) daughter's dR", 100,0,0.5);
    matchedDarkphotonMultiplicity_MuJet_= matchedDarkphotonDir.make<TH1F>("matchedDarkphotonMultiplicityMuJet", "Darkphoton(MuJet) Multiplicity",10,0,10);

    TFileDirectory jetInConeDir = fs_->mkdir("inCone_jet");
    jetInConeMultiplicity_ = jetInConeDir.make<TH1F>("jetInConeMultiplicity", "Jet in cone of darkphoton(not fully reconstructed) Multiplicity",5,0,5);
    jetInConeEnergy_       = jetInConeDir.make<TH1F>("jetInConeEnergy", "Jet Energy", 100,0,100);
    jetInConePt_           = jetInConeDir.make<TH1F>("jetInConePt", "Jet Pt", 100,0,100);
    jetInConeEt_           = jetInConeDir.make<TH1F>("jetInConeEt", "Jet Et", 100,0,100);
    jetInConeEta_          = jetInConeDir.make<TH1F>("jetInConeEta", "Jet Eta", 60,-3,3);
    jetInConePhi_          = jetInConeDir.make<TH1F>("jetInConePhi", "Jet Phi", 50,-3.1416,3.1416);
    jetInConeChargedHoverE_= jetInConeDir.make<TH1F>("jetInConeChargedHoverE", "Jet Charged H/E", 200,-1,1);
    jetInConeChargedMultiplicity_ = jetInConeDir.make<TH1F>("jetInConeChargedMultiplicity", "Jet Charged Multiplicity", 11,-1,10);
    jetInConeMuonEnergyFraction_  = jetInConeDir.make<TH1F>("jetInConeMuonEnergyFraction", "Jet Muon energy fraction", 200,-1,1);
    jetInConeMuonMultiplicity_    = jetInConeDir.make<TH1F>("jetInConeMuonMultiplicity", "Jet Muon Multiplicity",11,-1,10);
    jetInConeDR_                  = jetInConeDir.make<TH1F>("jetInConeDR", "dR between gen darkphoton and closest jet",100,0,0.5);
    jetInConeAloneMultiplicity_   = jetInConeDir.make<TH1F>("jetInConeAloneMultiplicity", "Jet in cone of darkphoton(no daughters reconstructed) Multiplicity",5,0,5);
    // eventTree_->Branch("numberOfElectrons", &electron_N, "numberOfElectrons/I");

}

// ------------ method called once each job just after ending the event loop  ------------
void 
sidm::MuonAnalysis::endJob() 
{
    muonMultiplicity_->GetXaxis()->SetTitle("Number of Muons");
    muonMultiplicity_->GetYaxis()->SetTitle("Event Number");
    
    muonEnergy_->GetXaxis()->SetTitle("Energy [GeV]");
    muonEnergy_->GetYaxis()->SetTitle("Event Number");

    muonPt_->GetXaxis()->SetTitle("pT [GeV]");
    muonPt_->GetYaxis()->SetTitle("Event Number");
    
    muonEt_->GetXaxis()->SetTitle("ET [GeV]");
    muonEt_->GetYaxis()->SetTitle("Event Number");

    muonEta_->GetXaxis()->SetTitle("eta");
    muonEta_->GetYaxis()->SetTitle("Event Number");
    
    muonPhi_->GetXaxis()->SetTitle("phi");
    muonPhi_->GetYaxis()->SetTitle("Event Number");
    
    muonVertex_->GetXaxis()->SetTitle("V_{z}");
    muonVertex_->GetYaxis()->SetTitle("V_{xy}");

    //------------------------------------------
    genMuonMultiplicity_->GetXaxis()->SetTitle("Number of Muons");
    genMuonMultiplicity_->GetYaxis()->SetTitle("Event Number");
    
    genMuonEnergy_->GetXaxis()->SetTitle("Energy [GeV]");
    genMuonEnergy_->GetYaxis()->SetTitle("Event Number");

    genMuonPt_->GetXaxis()->SetTitle("pT [GeV]");
    genMuonPt_->GetYaxis()->SetTitle("Event Number");

    genMuonEt_->GetXaxis()->SetTitle("ET [GeV]");
    genMuonEt_->GetYaxis()->SetTitle("Event Number");

    genMuonEta_->GetXaxis()->SetTitle("eta");
    genMuonEta_->GetYaxis()->SetTitle("Event Number");

    genMuonPhi_->GetXaxis()->SetTitle("phi");
    genMuonPhi_->GetYaxis()->SetTitle("Event Number");

    genMuonVertex_->GetXaxis()->SetTitle("V_{z}");
    genMuonVertex_->GetYaxis()->SetTitle("V_{xy}");

    //------------------------------------------
    matchedMuonMultiplicity_->GetXaxis()->SetTitle("Number of Muons");
    matchedMuonMultiplicity_->GetYaxis()->SetTitle("Event Number");
    
    matchedMuonEnergy_->GetXaxis()->SetTitle("Energy [GeV]");
    matchedMuonEnergy_->GetYaxis()->SetTitle("Event Number");

    matchedMuonPt_->GetXaxis()->SetTitle("pT [GeV]");
    matchedMuonPt_->GetYaxis()->SetTitle("Event Number");
    
    matchedMuonEt_->GetXaxis()->SetTitle("ET [GeV]");
    matchedMuonEt_->GetYaxis()->SetTitle("Event Number");

    matchedMuonEta_->GetXaxis()->SetTitle("eta");
    matchedMuonEta_->GetYaxis()->SetTitle("Event Number");
    
    matchedMuonPhi_->GetXaxis()->SetTitle("phi");
    matchedMuonPhi_->GetYaxis()->SetTitle("Event Number");
    
    matchedMuonVertex_->GetXaxis()->SetTitle("V_{z}");
    matchedMuonVertex_->GetYaxis()->SetTitle("V_{xy}");
    
    matchedMuonDR_->GetXaxis()->SetTitle("dR");
    matchedMuonDR_->GetYaxis()->SetTitle("Event Number");

    //------------------------------------------

    unmatchedMuonMultiplicity_->GetXaxis()->SetTitle("Number of Muons");
    unmatchedMuonMultiplicity_->GetYaxis()->SetTitle("Event Number");

    unmatchedMuonEnergy_->GetXaxis()->SetTitle("Energy [GeV]");
    unmatchedMuonEnergy_->GetYaxis()->SetTitle("Event Number");

    unmatchedMuonPt_->GetXaxis()->SetTitle("pT [GeV]");
    unmatchedMuonPt_->GetYaxis()->SetTitle("Event Number");

    unmatchedMuonEt_->GetXaxis()->SetTitle("ET [GeV]");
    unmatchedMuonEt_->GetYaxis()->SetTitle("Event Number");

    unmatchedMuonEta_->GetXaxis()->SetTitle("eta");
    unmatchedMuonEta_->GetYaxis()->SetTitle("Event Number");

    unmatchedMuonPhi_->GetXaxis()->SetTitle("phi");
    unmatchedMuonPhi_->GetYaxis()->SetTitle("Event Number");

    unmatchedMuonVertex_->GetXaxis()->SetTitle("V_{z}");
    unmatchedMuonVertex_->GetYaxis()->SetTitle("V_{xy}");

    //------------------------------------------

    matchedDarkphotonMultiplicity_MuMu_ ->GetXaxis()->SetTitle("Number of darkphotons");
    matchedDarkphotonMultiplicity_MuMu_ ->GetYaxis()->SetTitle("Event Number");

    matchedDarkphotonInvm_MuMu_         ->GetXaxis()->SetTitle("Invariant mass [GeV]");
    matchedDarkphotonInvm_MuMu_         ->GetYaxis()->SetTitle("Event Number");

    matchedDarkphotonDR_MuMu_           ->GetXaxis()->SetTitle("dR");
    matchedDarkphotonDR_MuMu_           ->GetXaxis()->SetTitle("Event Number");

    matchedDarkphotonMultiplicity_MuJet_ ->GetXaxis()->SetTitle("Number of darkphotons");
    matchedDarkphotonMultiplicity_MuJet_ ->GetYaxis()->SetTitle("Event Number");

    //------------------------------------------

    jetInConeMultiplicity_->GetXaxis()->SetTitle("Number of Jets");
    jetInConeMultiplicity_->GetYaxis()->SetTitle("Event number");

    jetInConeEnergy_      ->GetXaxis()->SetTitle("Energy [GeV]");
    jetInConeEnergy_      ->GetYaxis()->SetTitle("Event Number");

    jetInConePt_          ->GetXaxis()->SetTitle("pT [GeV]");
    jetInConePt_          ->GetYaxis()->SetTitle("Event Number");

    jetInConeEt_          ->GetXaxis()->SetTitle("ET [GeV]");
    jetInConeEt_          ->GetYaxis()->SetTitle("Event Number");
    
    jetInConeEta_          ->GetXaxis()->SetTitle("eta");
    jetInConeEta_          ->GetYaxis()->SetTitle("Event Number");
    
    jetInConePhi_          ->GetXaxis()->SetTitle("phi");
    jetInConePhi_          ->GetYaxis()->SetTitle("Event Number");

    jetInConeChargedHoverE_->GetXaxis()->SetTitle("Charged H/E");
    jetInConeChargedHoverE_->GetYaxis()->SetTitle("Event Number");

    jetInConeChargedMultiplicity_ ->GetXaxis()->SetTitle("Number of Charged particles");
    jetInConeChargedMultiplicity_ ->GetYaxis()->SetTitle("Event Number");

    jetInConeMuonEnergyFraction_  ->GetXaxis()->SetTitle("Muon energy fraction");
    jetInConeMuonEnergyFraction_  ->GetYaxis()->SetTitle("Event Number");

    jetInConeMuonMultiplicity_    ->GetXaxis()->SetTitle("Number of Muons");
    jetInConeMuonMultiplicity_    ->GetYaxis()->SetTitle("Event Number");

    jetInConeDR_                  ->GetXaxis()->SetTitle("dR");
    jetInConeDR_                  ->GetYaxis()->SetTitle("Event Number");

    jetInConeAloneMultiplicity_   ->GetXaxis()->SetTitle("Number of Jets");
    jetInConeAloneMultiplicity_   ->GetXaxis()->SetTitle("Event Number");
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
sidm::MuonAnalysis::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(sidm::MuonAnalysis);
