#ifndef SIDMBS_PHYSICSOBJECT_H
#define SIDMBS_PHYSICSOBJECT_H

#include <utility>
#include <cmath>
#include <vector>
#include "DataFormats/Common/interface/Ptr.h"
#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/Jet.h"

namespace sidm {

    class objBase{
    public:
        int _eventId;
        float _pt;
        float _eta;
        float _phi;
        float _energy;
        float _mass;
    };

    class objCompound : public objBase{
    public:
        float _invM;
        float _dEta;
        float _dPhi;
        float _dR;

    };
    
    class Jet : public objBase{
    public:
        Jet(){}
        Jet(const edm::Ptr<pat::Jet>& patj)
            : _charged_H_over_E( (patj->isPFJet() || patj->isJPTJet()) ? patj->chargedHadronEnergy()/patj->chargedEmEnergy() : -999 ),
              _H_over_E( patj->isCaloJet() ? patj->energyFractionHadronic()/patj->emEnergyFraction():-999 ),
              _chargedMultiplicity( (patj->isPFJet() || patj->isJPTJet()) ? patj->chargedMultiplicity():-999),
              _electronEnergyFraction( patj->isPFJet() ? patj->electronEnergyFraction() : -999 ),
              _electronMultiplicity( patj->isPFJet() ? patj->electronMultiplicity():-999),
              _num_elec_fullin( patj->isJPTJet() ?  patj->elecsInVertexInCalo().size():-999 ),
              _num_elec_curledin( patj->isJPTJet() ? patj->elecsOutVertexInCalo().size():-999 ),
              _num_elec_curledout( patj->isJPTJet() ? patj->elecsInVertexOutCalo().size():-999 )
        {
            _pt = patj->pt();
            _eta = patj->eta();
        }

        float _charged_H_over_E;  //< relative to uncorrected jet energy.
        float _H_over_E;
        int _chargedMultiplicity;
        float _electronEnergyFraction;
        int _electronMultiplicity;
        int _num_elec_fullin;
        int _num_elec_curledin;
        int _num_elec_curledout;
    };

    class Ep : public objBase{
    public:
        Ep(){}
        Ep(const Ep& ep_){}
        Ep(const edm::Ptr<reco::Candidate>& pate) {
            _pt     = pate->pt();
            _eta    = pate->eta();
            _phi    = pate->phi();
            _energy = pate->energy();
        }
        Ep(const reco::Candidate* gene) {
            _pt     = gene->pt();
            _eta    = gene->eta();
            _phi    = gene->phi();
            _energy = gene->energy();
        }
    };

    class Zp : public objCompound{
    public:
        Zp(){ matched = false; }
        Zp(const Ep& ine, const Ep& inp) : e(ine), p(inp) {
            matched = false;
            _dEta = std::abs(e._eta - p._eta);
            _dPhi = std::abs(e._phi - p._phi);
            _dR   = std::sqrt(_dEta*_dEta + _dPhi*_dPhi);
        }
        Zp(const std::pair<edm::Ptr<reco::Candidate>, edm::Ptr<reco::Candidate> >& q) {
            e = sidm::Ep(q.first);
            p = sidm::Ep(q.second);
            _dEta = std::abs(e._eta - p._eta);
            _dPhi = std::abs(e._phi - p._phi);
            _dR   = std::sqrt(_dEta*_dEta + _dPhi*_dPhi);
            _invM = (q.first->p4() + q.second->p4()).M();
            _eta  = e._eta + p._eta;
            _pt   = e._pt + p._pt;
        }
        Zp(const reco::Candidate* gene, const reco::Candidate* genp) {
            e = sidm::Ep(gene);
            p = sidm::Ep(genp);
            _dEta = std::abs(gene->eta() - genp->eta());
            _dPhi = std::abs(gene->phi() - genp->phi());
            _dR   = std::sqrt(_dEta*_dEta + _dPhi*_dPhi);
            _invM = (gene->p4() + genp->p4()).M();
            _eta  = e._eta + p._eta;
            _pt   = e._pt + p._pt;
            _dv_x = gene->vx();
            _dv_y = gene->vy();
            _dv_z = gene->vz();
        }
        
        Ep e, p;
        bool matched;
        float _dv_x;
        float _dv_y;
        float _dv_z;
        std::pair<float, float> dRVal (const std::pair<edm::Ptr<reco::Candidate>, edm::Ptr<reco::Candidate> >& q) const {
            std::pair<float, float> tmp(0., 0.);
            tmp.first = std::sqrt( (e._eta-q.first->eta())*(e._eta-q.first->eta()) + (e._phi-q.first->phi())*(e._phi-q.first->phi()) );
            tmp.second = std::sqrt( (p._eta-q.second->eta())*(p._eta-q.second->eta()) + (p._phi-q.second->phi())*(p._phi-q.second->phi()) );
            return tmp;
        }
        std::pair<float, float> dRVal (const edm::Ptr<pat::Jet>& j) const {
            std::pair<float, float> tmp(0., 0.);
            tmp.first  = std::sqrt( (e._eta-j->eta())*(e._eta-j->eta()) + (e._phi-j->phi())*(e._phi-j->phi()) );
            tmp.second = std::sqrt( (p._eta-j->eta())*(p._eta-j->eta()) + (p._phi-j->phi())*(p._phi-j->phi()) );
            return tmp;
        }
    };

    class Ps : public objCompound{
    public:
        Zp zp_0;
        Zp zp_1;
    };

} // namespace sidm

#endif
