// MEMBER 1 | tiny test helper (no external libraries needed)
#pragma once
#include <iostream>
inline int g_pass = 0, g_fail = 0;
#define CHECK(cond)                                                              \
    do {                                                                         \
        if (cond) { ++g_pass; }                                                  \
        else { ++g_fail; std::cerr << "  FAIL " << __FILE__ << ":" << __LINE__   \
                                   << "  " << #cond << "\n"; }                   \
    } while (0)
#define SUMMARY()                                                                \
    do {                                                                         \
        std::cout << g_pass << " passed, " << g_fail << " failed\n";             \
        return g_fail ? 1 : 0;                                                   \
    } while (0)
