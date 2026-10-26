// =====================================================================
// MEMBER 1  |  dna_metrics.hpp
// DNA helpers + distance functions (the "metric" the VP tree relies on)
// =====================================================================
#pragma once
#include <algorithm>
#include <cctype>
#include <string>
#include <vector>

namespace dna {

// A distance function: must be a true metric (d>=0, d(a,a)=0,
// symmetric, triangle inequality) or VP-tree pruning is not valid.
using Metric = int (*)(const std::string&, const std::string&);

// One search result: index of the read in the database + its distance
struct Match {
    int index;
    int dist;
};

inline bool matchLess(const Match& a, const Match& b) {
    return a.dist < b.dist || (a.dist == b.dist && a.index < b.index);
}

// Global counter of how many distance computations happened.
// Used in the demo to PROVE the VP tree prunes work vs brute force.
inline long long& distanceCalls() {
    static long long c = 0;
    return c;
}
inline void resetDistanceCalls() { distanceCalls() = 0; }

// Uppercase and remove whitespace
inline std::string normalize(const std::string& s) {
    std::string out;
    out.reserve(s.size());
    for (char c : s)
        if (!std::isspace(static_cast<unsigned char>(c)))
            out.push_back(static_cast<char>(std::toupper(static_cast<unsigned char>(c))));
    return out;
}

// Only A, C, G, T, N allowed
inline bool isValid(const std::string& s) {
    for (char c : s)
        if (c != 'A' && c != 'C' && c != 'G' && c != 'T' && c != 'N') return false;
    return true;
}

inline std::string reverseComplement(const std::string& s) {
    std::string r(s.rbegin(), s.rend());
    for (char& c : r) {
        switch (c) {
            case 'A': c = 'T'; break;
            case 'T': c = 'A'; break;
            case 'C': c = 'G'; break;
            case 'G': c = 'C'; break;
            default: break;  // N stays N
        }
    }
    return r;
}

// Hamming distance, O(min(m,n)). Length difference counts as mismatches
// (equivalent to padding the shorter read with a unique symbol) so it
// stays a valid metric even for reads of unequal length.
inline int hamming(const std::string& a, const std::string& b) {
    ++distanceCalls();
    size_t n = std::min(a.size(), b.size());
    int d = static_cast<int>(std::max(a.size(), b.size()) - n);
    for (size_t i = 0; i < n; ++i)
        if (a[i] != b[i]) ++d;
    return d;
}

// Levenshtein (edit) distance: min insertions+deletions+substitutions.
// Time O(m*n), space O(min(m,n)) using two rolling rows.
inline int levenshtein(const std::string& a, const std::string& b) {
    ++distanceCalls();
    const std::string& longer = a.size() >= b.size() ? a : b;
    const std::string& shorter = a.size() >= b.size() ? b : a;
    std::vector<int> prev(shorter.size() + 1), cur(shorter.size() + 1);
    for (size_t j = 0; j <= shorter.size(); ++j) prev[j] = static_cast<int>(j);
    for (size_t i = 1; i <= longer.size(); ++i) {
        cur[0] = static_cast<int>(i);
        for (size_t j = 1; j <= shorter.size(); ++j) {
            int cost = (longer[i - 1] == shorter[j - 1]) ? 0 : 1;
            cur[j] = std::min({prev[j] + 1, cur[j - 1] + 1, prev[j - 1] + cost});
        }
        std::swap(prev, cur);
    }
    return prev[shorter.size()];
}

inline const char* metricName(Metric m) {
    return m == &hamming ? "Hamming" : "Levenshtein (edit distance)";
}

}  // namespace dna
