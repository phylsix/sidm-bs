#ifndef SIDMBS_UTILITIES_H
#define SIDMBS_UTILITIES_H

#include <vector>
#include <utility>
#include <algorithm>

#include "DataFormats/Common/interface/Ptr.h"
#include "DataFormats/Candidate/interface/Candidate.h"

#include "sidm-bs/jetComponent/interface/physicsObject.h"

namespace sidm {

    class test {
    public:
        test(){};
        test(int a) : _a(a) {};
        virtual ~test() {};
        void get();
    private:
        int _a;

    };
    
    template<class T, class U=T> class pairvec {
    using pv = std::vector<std::pair<edm::Ptr<T>, edm::Ptr<U> > >;
    public:
        pairvec(const std::vector<edm::Ptr<T> >& a,
                const std::vector<edm::Ptr<U> >& b) :
            _a(a), _b(b) {};


        pairvec(const pv& vp, bool unique=false) {
            for (const auto& p : vp) {
                if (!unique){
                    _a.push_back(p.first);
                    _b.push_back(p.second);
                } else {
                    if (std::find(cbegin(_a), cend(_a), p.first) == cend(_a) &&
                        std::find(cbegin(_b), cend(_b), p.second)== cend(_b) ) {
                        _a.push_back(p.first);
                        _b.push_back(p.second);
                    }
                }
            }
        }
        

        virtual ~pairvec() {};


        // return a vector whose size = _a.size()*_b.size()
        pv get() const {
            pv tmp{};
            for (const auto& p : _a) {
                for (const auto& q : _b) {
                    tmp.emplace_back(std::make_pair(p, q));
                }
            }
            return tmp;
        }


        // zip pairs, requires _a.size() == _b.size().
        // e.g.[(_a[0], _b[0]),(_a[1], _b[1]),...]
        pv get_zip() const {
            assert(_a.size() == _b.size());
            pv tmp{};
            typename std::vector<edm::Ptr<T> >::const_iterator ai;
            typename std::vector<edm::Ptr<U> >::const_iterator bi;
            for ( ai=_a.begin(), bi=_b.begin();
                    ai!=_a.end(), bi!=_b.end();
                    ++ai, ++bi) {
                tmp.emplace_back(std::make_pair(*ai, *bi));
            }
            return tmp;
        }

    private:
        std::vector<edm::Ptr<T> > _a;
        std::vector<edm::Ptr<U> > _b;
    };

    void match_patPair_with_zps (const std::pair<edm::Ptr<reco::Candidate>, edm::Ptr<reco::Candidate> >& undertest,
                                 std::vector<sidm::Zp>& zps,
                                 double deltaR_upsideband);

    template<class T>
    void remove_from_collection (std::vector<T>* collection, const T& val) {
        auto to_be_removed = std::find(collection->begin(), collection->end(), val);
        if (to_be_removed != collection->end())
            collection->erase(to_be_removed);
    }
}

#endif
