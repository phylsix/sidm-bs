#ifndef SIDMBS_PHYSICSOBJECT_H
#define SIDMBS_PHYSICSOBJECT_H

#include <utility>
#include <cmath>
#include "DataFormats/Common/interface/Ptr.h"
#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/PatCandidates/interface/Electron.h"

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
    };

    class Ep : public objBase{
    public:
        Ep(){}
        Ep(const Ep& ep_){}
        Ep(edm::Ptr<pat::Electron> pate) {
            _pt     = pate->pt();
            _eta    = pate->eta();
            _phi    = pate->phi();
            _energy = pate->energy();
        }
    };

    class Zp : public objCompound{
    public:
        Zp(){ matched = false; }
        Zp(Ep& ine, Ep& inp) : e(ine), p(inp) {
            matched = false;
            _dEta = std::abs(e._eta - p._eta);
            _dPhi = std::abs(e._phi - p._phi);
            _dR   = std::sqrt(_dEta*_dEta + _dPhi*_dPhi);
        }
        Zp(std::pair<edm::Ptr<pat::Electron>, edm::Ptr<pat::Electron> > q) {
            e = sidm::Ep(q.first);
            p = sidm::Ep(q.second);
            _dEta = std::abs(e._eta - p._eta);
            _dPhi = std::abs(e._phi - p._phi);
            _dR   = std::sqrt(_dEta*_dEta + _dPhi*_dPhi);
        }
        
        Ep e, p;
        bool matched;
        float _dv_x;
        float _dv_y;
        float _dv_z;
        std::pair<float, float> dRVal (std::pair<edm::Ptr<reco::Candidate>, edm::Ptr<reco::Candidate> > q) const {
            std::pair<float, float> tmp(0., 0.);
            tmp.first = std::sqrt( (e._eta-q.first->eta())*(e._eta-q.first->eta()) + (e._phi-q.first->phi())*(e._phi-q.first->phi()) );
            tmp.second = std::sqrt( (p._eta-q.second->eta())*(p._eta-q.second->eta()) + (p._phi-q.second->phi())*(p._phi-q.second->phi()) );
            return tmp;
        }
    };

    class Ps : public objCompound{
    public:
        Zp zp_0;
        Zp zp_1;
    };
}

#endif
