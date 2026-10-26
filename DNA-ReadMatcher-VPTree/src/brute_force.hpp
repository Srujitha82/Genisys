// =====================================================================
// MEMBER 1  |  brute_force.hpp
// Linear-scan baseline. Used (a) to verify the VP tree is correct and
// (b) to compare speed / number of distance computations.
// =====================================================================
#pragma once
#include <algorithm>
#include <string>
#include <vector>
#include "dna_metrics.hpp"

namespace dna {

inline std::vector<Match> bruteRange(const std::vector<std::string>& reads,
                                     const std::string& q, int radius, Metric m) {
    std::vector<Match> out;
    for (size_t i = 0; i < reads.size(); ++i) {
        int d = m(q, reads[i]);
        if (d <= radius) out.push_back({static_cast<int>(i), d});
    }
    std::sort(out.begin(), out.end(), matchLess);
    return out;
}

inline std::vector<Match> bruteKNN(const std::vector<std::string>& reads,
                                   const std::string& q, int k, Metric m) {
    std::vector<Match> all;
    for (size_t i = 0; i < reads.size(); ++i)
        all.push_back({static_cast<int>(i), m(q, reads[i])});
    std::sort(all.begin(), all.end(), matchLess);
    if (k < 0) k = 0;
    if (static_cast<size_t>(k) < all.size()) all.resize(k);
    return all;
}

}  // namespace dna
