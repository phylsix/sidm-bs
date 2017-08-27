#include "../interface/utilities.h"

#include <iostream>
#include <algorithm>

void
sidm::test::get() {
    std::cout << "test: " << _a << std::endl;
}

/*--
template<class T>
sidm::pairvec<T>::pairvec(std::vector<std::pair<edm::Ptr<T>, edm::Ptr<T> > >& vp, bool unique) {
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

template<class T>
std::vector<std::pair<edm::Ptr<T>, edm::Ptr<T> > >
sidm::pairvec<T>::get() {
    //std::vector<std::pair<edm::Ptr<T>, edm::Ptr<T> > > tmp{};
    for (const auto& p : _a) {
        for (const auto& q : _b) {
            n.emplace_back(std::make_pair(p, q));
        }
    }
    return n;
}

template<class T>
std::vector<std::pair<edm::Ptr<T>, edm::Ptr<T> > >
sidm::pairvec<T>::get_zip() const {

    assert(_a.size() == _b.size());
    //std::vector<std::pair<edm::Ptr<T>, edm::Ptr<T> > > tmp{};
    typename std::vector<edm::Ptr<T> >::const_iterator ai, bi;
    for ( ai=_a.begin(), bi=_b.begin();
          ai!=_a.end(), bi!=_b.end();
          ++ai, ++bi) {
        z.emplace_back(std::make_pair(*ai, *bi));
    }

    return z;

}
--*/
