#include "../interface/utilities.h"

#include <iostream>
#include <algorithm>

void
sidm::test::get() {
    std::cout << "test: " << _a << std::endl;
}

/**
 * Test a pair of pat::Electron against two gen Zps with dR limit.
 * set .matched property of sidm::Zp if matched.
 *
 * @param [const std::pair<edm::Ptr<reco::Candidate>, edm::Ptr<reco::Candidate> >&] under_test pair of electrons to be tested
 * @param [std::vector<sidm::Zp>&] zps vector(2) gen zps tested against
 * @param [double] deltaR_upsideband limit of deltaR
 */
void
sidm::match_patPair_with_zps (const std::pair<edm::Ptr<reco::Candidate>, edm::Ptr<reco::Candidate> >& under_test,
                              std::vector<sidm::Zp>& zps,
                              double deltaR_upsideband) {
    const int matchedNum = std::count_if(std::cbegin(zps), std::cend(zps),
            [](const sidm::Zp& p){return p.matched;});

    if ( matchedNum == 0 ) {
        std::pair<float, float> dr_tmp0(zps[0].dRVal(under_test));
        std::pair<float, float> dr_tmp1(zps[1].dRVal(under_test));
        float average_dr0 = (dr_tmp0.first + dr_tmp0.second)/2;
        float average_dr1 = (dr_tmp1.first + dr_tmp1.second)/2;
        if (average_dr0<=deltaR_upsideband && average_dr1>deltaR_upsideband)
            zps[0].matched = true;
        if (average_dr1<=deltaR_upsideband && average_dr0>deltaR_upsideband)
            zps[1].matched = true;
        if (average_dr0<=deltaR_upsideband && average_dr1<=deltaR_upsideband) {
            if (average_dr0<=average_dr1) 
                zps[0].matched = true;
            else
                zps[1].matched = true;
        }
    
    } else if ( matchedNum == 1 ) {
        for (auto& p : zps) {
            if (p.matched) continue;
            std::pair<float, float> dr_tmp(p.dRVal(under_test));
            if ((dr_tmp.first+dr_tmp.second)/2 < deltaR_upsideband)
                p.matched = true;
        }
    }
}

/*--
void
sidm::match_jet_with_zps (const edm::Ptr<pat::Jet>& under_test,
                          std::vector<sidm::Zp>& zps,
                          double deltaR_upsideband) {
    std::vector<float> average_deltaR{};
    for (auto& p : zps) {
        if (p.matched) continue;
        std::pair<float, float> delta_tmp(p.dRVal(under_test));
        if (delta_tmp.first <= deltaR_upsideband && delta_tmp.second <= deltaR_upsideband)
            p.matched = true;
        else
            p.matched = false;
        average_deltaR.push_back((delta_tmp.first+delta_tmp.second)/2);
    }

    int matchedNum = std::count_if(std::cbegin(zps), std::cend(zps),
            [](const sidm::Zp& p){return p.matched;} );

    if (matched==2 && delta_tmp.size()==2) {
        // Whose average dR is smaller got matched,
        // the other one reamains as unmatched.
        if (average_deltaR[0] <= average_deltaR[1]) 
            zps[1].matched = false;
        else
            zps[0].matched = false;
    }
}
--*/
