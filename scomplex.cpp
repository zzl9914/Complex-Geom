/**
 * @file scomplex.cpp
 * @brief Non-template implementations: round_left, medium_mod,
 *        shift_left, shift_right.
 *
 * Templates that take a user functor (differential, invert, integral)
 * and the rest of the class/function templates stay in scomplex.hpp.
 * Callers must compile and link this file.
 */

#include "scomplex.hpp"

namespace std {
namespace complex {
namespace experimental {
namespace zhangzl {

double round_left(double a) { return a - std::round(a); }

double& medium_mod(double& a, const double b) {
    a /= b;
    a = round_left(a);
    a *= b;
    return a;
}

double shift_left(double a, uint64_t b) noexcept {
    if (!a) return 0;
    union {
        double a;
        uint64_t ap;
    } _a;
    _a.a = a;
    _a.ap += b << 52;
    return _a.a;
}

double shift_right(double a, uint64_t b) noexcept {
    if (!a) return 0;
    union {
        double a;
        uint64_t ap;
    } _a;
    _a.a = a;
    _a.ap -= b << 52;
    return _a.a;
}

}
}
}
}
