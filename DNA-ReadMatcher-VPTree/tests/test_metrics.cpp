// MEMBER 1 | unit tests for distance functions and DNA helpers
#include "../src/dna_metrics.hpp"
#include "../src/brute_force.hpp"
#include "check.hpp"
using namespace dna;

int main() {
    // known values
    CHECK(levenshtein("KITTEN", "SITTING") == 3);
    CHECK(levenshtein("ACGT", "ACGT") == 0);
    CHECK(levenshtein("", "ACGT") == 4);
    CHECK(levenshtein("ACGT", "") == 4);
    CHECK(levenshtein("", "") == 0);
    CHECK(levenshtein("ACGT", "AGT") == 1);   // one deletion
    CHECK(levenshtein("ACGT", "AGGT") == 1);  // one substitution
    CHECK(hamming("ACGT", "AGGT") == 1);
    CHECK(hamming("AAAA", "TTTT") == 4);
    CHECK(hamming("ACGT", "AC") == 2);        // length difference counted
    CHECK(hamming("", "") == 0);

    // metric axioms on a sample set
    std::vector<std::string> s = {"", "A", "ACGT", "AGGT", "TTTT", "ACGTACGT", "GATTACA", "N"};
    Metric ms[2] = {&levenshtein, &hamming};
    for (Metric m : ms)
        for (auto& a : s)
            for (auto& b : s) {
                CHECK(m(a, b) == m(b, a));               // symmetry
                CHECK((m(a, b) == 0) == (a == b));       // identity
                for (auto& c : s) CHECK(m(a, c) <= m(a, b) + m(b, c));  // triangle
            }

    // helpers
    CHECK(normalize(" ac gt\n") == "ACGT");
    CHECK(isValid("ACGTN"));
    CHECK(!isValid("ACGU"));
    CHECK(reverseComplement("AACG") == "CGTT");

    // brute force baseline
    std::vector<std::string> db = {"AAAA", "AAAT", "TTTT"};
    auto r = bruteRange(db, "AAAA", 1, &levenshtein);
    CHECK(r.size() == 2 && r[0].index == 0 && r[1].index == 1);
    auto k = bruteKNN(db, "TTTT", 1, &levenshtein);
    CHECK(k.size() == 1 && k[0].index == 2);
    CHECK(bruteKNN(db, "A", 0, &levenshtein).empty());
    CHECK(bruteKNN(db, "A", 99, &levenshtein).size() == 3);

    SUMMARY();
}
