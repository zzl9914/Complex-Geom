/**
 * @file scomplex.c
 * @brief Implementation of the C API declared in scomplex.h.
 *
 * Internal Newton / AGM / hyperoperation helpers are file-static.
 * sc_last_error is a single process-wide object shared by every
 * translation unit that links this file.
 */

#include "scomplex.h"

#include <ctype.h>
#include <float.h>
#include <limits.h>
#include <math.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>

sc_error_t sc_last_error = SC_OK;

const char* sc_error_string(sc_error_t e) {
    switch (e) {
        case SC_OK:              return "Success";
        case SC_ERR_DIV_BY_ZERO: return "Division by zero in complex number";
        case SC_ERR_OVERFLOW:    return "Overflow in complex number operation";
        case SC_ERR_DOMAIN:      return "Domain error in complex operation";
        case SC_ERR_CONVERGENCE: return "Failed to converge in complex operation";
        case SC_ERR_INVALID:     return "Invalid argument";
        default:                 return "Unknown error";
    }
}

void sc_set_error(sc_error_t e) { sc_last_error = e; }

sc_error_t sc_get_error(void) {
    sc_error_t e = sc_last_error;
    sc_last_error = SC_OK;
    return e;
}

sc_error_t sc_peek_error(void) { return sc_last_error; }

static void sc_check_overflow(double v) {
    if (!isfinite(v)) sc_set_error(SC_ERR_OVERFLOW);
    if (v < 0) v = -v;
    if (v > DBL_MAX / 2.0) sc_set_error(SC_ERR_OVERFLOW);
}

static void sc_check_div_zero(double d) {
    if (d < 0) d = -d;
    if (d < DBL_EPSILON) sc_set_error(SC_ERR_DIV_BY_ZERO);
}

typedef union { float value;  int32_t digits; } sc_ieee754_f;

typedef union { double value; uint64_t digits; } sc_ieee754_d;

float sc_flnot(float a) { sc_ieee754_f u; u.value = a; u.digits = ~u.digits; return u.value; }

float sc_flor(float a, float b) {
    sc_ieee754_f u, v; u.value = a; v.value = b; u.digits |= v.digits; return u.value;
}

float sc_fland(float a, float b) { return sc_flnot(sc_flor(sc_flnot(a), sc_flnot(b))); }

float sc_flxor(float a, float b) { return sc_fland(sc_flor(a, b), sc_flnot(sc_fland(a, b))); }

double sc_dblnot(double a) { sc_ieee754_d u; u.value = a; u.digits = ~u.digits; return u.value; }

double sc_dblor(double a, double b) {
    sc_ieee754_d u, v, ans;
    u.value = a; v.value = b;
    ans.digits = u.digits | v.digits;
    if (!(ans.digits & ((1ULL << 52) - 1))) return 0.0;
    return ans.value;
}

double sc_dbland(double a, double b) { return sc_dblnot(sc_dblor(sc_dblnot(a), sc_dblnot(b))); }

double sc_dblxor(double a, double b) { return sc_dbland(sc_dblor(a, b), sc_dblnot(sc_dbland(a, b))); }

double sc_shift_left(double a, uint64_t b) {
    if (a == 0.0) return 0.0;
    sc_ieee754_d u;
    u.value = a;
    u.digits += b << 52;
    return u.value;
}

double sc_shift_right(double a, uint64_t b) {
    if (a == 0.0) return 0.0;
    sc_ieee754_d u;
    u.value = a;
    u.digits -= b << 52;
    return u.value;
}

double sc_round_left(double a) { return a - round(a); }

double sc_pi(void) { return 3.14159265358979323846; }

double sc_2pi(void) { return 6.28318530717958647692; }

double sc_pi_2(void) { return 1.57079632679489661923; }

double sc_pi_3(void) { return 1.04719755119659774615; }

double sc_e(void) { return 2.71828182845904523536; }

double sc_sqrt2(void) { return 1.41421356237309504880; }

double sc_sqrt1_2(void) { return 0.70710678118654752440; }

double sc_cbrt2(void) { return 1.25992104989487316476; }

double sc_cbrt1_2(void) { return 0.79370052598409973737; }

sc_complex sc_make(double real, double imag) {
    sc_check_overflow(real);
    sc_check_overflow(imag);
    sc_complex c;
    c.real = real;
    c.imag = imag;
    return c;
}

sc_complex sc_make_real(double real) { return sc_make(real, 0.0); }

sc_complex sc_zero(void) { return sc_make(0.0, 0.0); }

sc_complex sc_one(void) { return sc_make(1.0, 0.0); }

sc_complex sc_i(void) { return sc_make(0.0, 1.0); }

sc_complex sc_exp_i_pi_1_6(void) {
    return sc_make(0.5, 0.86602540378443864676);
}

sc_complex sc_exp_i_pi_1_9(void) {
    return sc_make(0.93969262078590838405, 0.34202014332566873304);
}

sc_complex sc_from_string(const char* str) {
    char buf[256];
    size_t n = 0;
    for (const char* s = str; *s && n + 1 < sizeof(buf); ++s) {
        if (!isspace((unsigned char)*s)) buf[n++] = *s;
    }
    buf[n] = 0;
    if (n == 0) { sc_set_error(SC_ERR_INVALID); return sc_zero(); }

    if (strcmp(buf, "i") == 0 || strcmp(buf, "+i") == 0) return sc_i();
    if (strcmp(buf, "-i") == 0) return sc_make(0.0, -1.0);

    char* end = NULL;
    double real = strtod(buf, &end);
    if (end != buf && *end == 0) return sc_make(real, 0.0);
    if (end == buf) { sc_set_error(SC_ERR_INVALID); return sc_zero(); }

    char sign = *end;
    if (sign != '+' && sign != '-') { sc_set_error(SC_ERR_INVALID); return sc_zero(); }
    const char* q = end + 1;
    double imag = 1.0;
    if (*q != 'i' && *q != 'I') {
        char* end2 = NULL;
        imag = strtod(q, &end2);
        if (end2 == q) { sc_set_error(SC_ERR_INVALID); return sc_zero(); }
        q = end2;
    }
    if (*q != 'i' && *q != 'I') { sc_set_error(SC_ERR_INVALID); return sc_zero(); }
    if (q[1] != 0) { sc_set_error(SC_ERR_INVALID); return sc_zero(); }
    if (sign == '-') imag = -imag;
    return sc_make(real, imag);
}

double sc_real(sc_complex a) { return a.real; }

double sc_imag(sc_complex a) { return a.imag; }

bool sc_equal(sc_complex a, sc_complex b) { return a.real == b.real && a.imag == b.imag; }

bool sc_not_equal(sc_complex a, sc_complex b) { return !sc_equal(a, b); }

bool sc_less(sc_complex a, sc_complex b) {
    return (a.real == b.real) ? (a.imag < b.imag) : (a.real < b.real);
}

bool sc_greater(sc_complex a, sc_complex b) { return !sc_less(a, b) && !sc_equal(a, b); }

bool sc_less_equal(sc_complex a, sc_complex b) { return sc_less(a, b) || sc_equal(a, b); }

bool sc_greater_equal(sc_complex a, sc_complex b) { return !sc_less(a, b); }

bool sc_is_zero(sc_complex a) { return a.real == 0.0 && a.imag == 0.0; }

bool sc_nonzero(sc_complex a) { return !sc_is_zero(a); }

bool sc_is_finite(sc_complex a) { return isfinite(a.real) && isfinite(a.imag); }

bool sc_is_nan(sc_complex a) { return isnan(a.real) || isnan(a.imag); }

sc_complex sc_add(sc_complex a, sc_complex b) {
    return sc_make(a.real + b.real, a.imag + b.imag);
}

sc_complex sc_sub(sc_complex a, sc_complex b) {
    return sc_make(a.real - b.real, a.imag - b.imag);
}

sc_complex sc_mul(sc_complex a, sc_complex b, double system) {
    return sc_make(a.real * b.real + system * a.imag * b.imag,
                   a.imag * b.real + a.real * b.imag);
}

sc_complex sc_div(sc_complex a, sc_complex b, double system) {
    double denom = b.real * b.real - system * b.imag * b.imag;
    sc_check_div_zero(denom);
    if (sc_peek_error() == SC_ERR_DIV_BY_ZERO) return sc_zero();
    return sc_make((a.real * b.real - system * a.imag * b.imag) / denom,
                   (a.imag * b.real - a.real * b.imag) / denom);
}

sc_complex sc_mod(sc_complex a, sc_complex b, double system) {
    sc_complex q = sc_div(a, b, system);
    q.real = sc_round_left(q.real);
    q.imag = sc_round_left(q.imag);
    return sc_mul(q, b, system);
}

sc_complex sc_neg(sc_complex a) { return sc_make(-a.real, -a.imag); }

sc_complex sc_conj(sc_complex a) { return sc_make(a.real, -a.imag); }

sc_complex sc_add_real(sc_complex a, double b) { return sc_make(a.real + b, a.imag); }

sc_complex sc_sub_real(sc_complex a, double b) { return sc_make(a.real - b, a.imag); }

sc_complex sc_real_sub(double a, sc_complex b) { return sc_make(a - b.real, -b.imag); }

sc_complex sc_mul_real(sc_complex a, double b) { return sc_make(a.real * b, a.imag * b); }

sc_complex sc_div_real(sc_complex a, double b) {
    sc_check_div_zero(b);
    if (sc_peek_error() == SC_ERR_DIV_BY_ZERO) return sc_zero();
    return sc_make(a.real / b, a.imag / b);
}

sc_complex sc_real_div(double a, sc_complex b, double system) {
    return sc_div(sc_make_real(a), b, system);
}

sc_complex sc_mod_real(sc_complex a, double b) {
    sc_check_div_zero(b);
    if (sc_peek_error() == SC_ERR_DIV_BY_ZERO) return sc_zero();
    return sc_make(sc_round_left(a.real / b) * b, sc_round_left(a.imag / b) * b);
}

sc_complex sc_real_mod(double a, sc_complex b, double system) {
    return sc_mod(sc_make_real(a), b, system);
}

sc_complex sc_shl(sc_complex a, uint64_t b) {
    return sc_make(sc_shift_left(a.real, b), sc_shift_left(a.imag, b));
}

sc_complex sc_shr(sc_complex a, uint64_t b) {
    return sc_make(sc_shift_right(a.real, b), sc_shift_right(a.imag, b));
}

double sc_norm(sc_complex a, double system) {
    return a.real * a.real - system * a.imag * a.imag;
}

double sc_norm_norm(sc_complex a) {
    return a.real * a.real + a.imag * a.imag;
}

sc_complex sc_abs(sc_complex a) {
    return sc_make(fabs(a.real), fabs(a.imag));
}

sc_complex sc_fabs(sc_complex a) { return sc_abs(a); }

double sc_abs2(sc_complex a) { return sc_norm_norm(a); }

double sc_modulus(sc_complex a) { return sqrt(sc_norm_norm(a)); }

static int sc_failed(void) { return sc_peek_error() != SC_OK; }

sc_complex sc_pow_int(sc_complex base, int exponent, double system) {
    if (exponent == 0) return sc_one();
    if (sc_is_zero(base)) return sc_zero();
    /* Iterative. Negating INT_MIN does not change the value, so the old
       recursive `1 / pow(base, -exponent)` never returned. */
    int neg = exponent < 0;
    unsigned uexp = neg
        ? (exponent == INT_MIN ? 2147483648u : (unsigned)(-exponent))
        : (unsigned)exponent;
    sc_complex result = sc_one();
    sc_complex factor = base;
    while (uexp) {
        if (uexp & 1u) result = sc_mul(result, factor, system);
        uexp >>= 1;
        if (uexp) factor = sc_mul(factor, factor, system);
    }
    if (neg) return sc_div(sc_one(), result, system);
    return result;
}

static sc_complex sc_sqrt_init(sc_complex a, sc_complex b, double system, int max_iterations) {
    sc_complex bp = b;
    for (int i = 0; i < max_iterations; ++i) {
        sc_error_t saved = sc_peek_error();
        sc_complex nx = sc_div_real(sc_add(bp, sc_div(a, bp, system)), 2.0);
        if (sc_peek_error() != SC_OK && sc_peek_error() != saved) {
            sc_set_error(SC_ERR_CONVERGENCE);
            return b;
        }
        sc_complex d = sc_abs(sc_sub(nx, bp));
        b = nx;
        if (d.real < DBL_EPSILON && d.imag < DBL_EPSILON) break;
        bp = b;
    }
    return b;
}

static sc_complex sc_sqrt_iter(sc_complex a, double system, int max_iterations) {
    if (a.imag == 0.0 && a.real >= 0.0) return sc_make(sqrt(a.real), 0.0);
    if (sc_norm(a, system) < 0.0) { sc_set_error(SC_ERR_DOMAIN); return sc_zero(); }
    double L = sc_norm_norm(a);
    if (L > 2.0) return sc_mul_real(sc_sqrt_iter(sc_shr(a, 1), system, max_iterations), sc_sqrt2());
    if (L < 0.5) return sc_mul_real(sc_sqrt_iter(sc_shl(a, 1), system, max_iterations), sc_sqrt1_2());
    if (a.real > 0.0) return sc_sqrt_init(a, sc_one(), system, max_iterations);
    if (a.imag > 0.0) return sc_sqrt_init(a, sc_exp_i_pi_1_6(), system, max_iterations);
    return sc_sqrt_init(a, sc_conj(sc_exp_i_pi_1_6()), system, max_iterations);
}

sc_complex sc_sqrt(sc_complex a, double system) { return sc_sqrt_iter(a, system, 100); }

static sc_complex sc_cbrt_init(sc_complex a, sc_complex b, double system, int max_iterations) {
    sc_complex bp = b;
    for (int i = 0; i < max_iterations; ++i) {
        sc_error_t saved = sc_peek_error();
        sc_complex nx = sc_div_real(sc_add(sc_mul_real(bp, 2.0),
                                           sc_div(a, sc_mul(bp, bp, system), system)), 3.0);
        if (sc_peek_error() != SC_OK && sc_peek_error() != saved) {
            sc_set_error(SC_ERR_CONVERGENCE);
            return b;
        }
        sc_complex d = sc_abs(sc_sub(nx, bp));
        b = nx;
        if (d.real < DBL_EPSILON && d.imag < DBL_EPSILON) break;
        bp = b;
    }
    return b;
}

static sc_complex sc_cbrt_iter(sc_complex a, double system, int max_iterations) {
    if (a.imag == 0.0) return sc_make(cbrt(a.real), 0.0);
    double L = sc_norm_norm(a);
    if (L > 2.0) return sc_mul_real(sc_cbrt_iter(sc_shr(a, 1), system, max_iterations), sc_cbrt2());
    if (L < 0.5) return sc_mul_real(sc_cbrt_iter(sc_shl(a, 1), system, max_iterations), sc_cbrt1_2());
    if (a.real > 0.0) return sc_cbrt_init(a, sc_one(), system, max_iterations);
    if (a.imag > 0.0) return sc_cbrt_init(a, sc_exp_i_pi_1_9(), system, max_iterations);
    return sc_cbrt_init(a, sc_conj(sc_exp_i_pi_1_9()), system, max_iterations);
}

sc_complex sc_cbrt(sc_complex a, double system) { return sc_cbrt_iter(a, system, 100); }

static sc_complex sc_root_iter(sc_complex a, int n, double system, int max_iterations) {
    if (n == 1) return a;
    if (n == 2) return sc_sqrt_iter(a, system, max_iterations);
    if (n == 3) return sc_cbrt_iter(a, system, max_iterations);
    if (a.imag == 0.0) return sc_make(pow(a.real, 1.0 / n), 0.0);
    double L = sc_norm_norm(a);
    if (L > 2.0) return sc_mul_real(sc_root_iter(sc_shr(a, 1), n, system, max_iterations), pow(2.0, 1.0 / n));
    if (L < 0.5) return sc_mul_real(sc_root_iter(sc_shl(a, 1), n, system, max_iterations), pow(0.5, 1.0 / n));
    sc_complex init = sc_make(cos(sc_pi_3() / n), sin(sc_pi_3() / n));
    if (a.real > 0.0) init = sc_one();
    else if (a.imag <= 0.0) init = sc_conj(init);
    sc_complex b = init, bp = init;
    for (int i = 0; i < max_iterations; ++i) {
        sc_error_t saved = sc_peek_error();
        sc_complex nx = sc_div_real(sc_add(sc_mul_real(bp, (double)(n - 1)),
                                           sc_div(a, sc_pow_int(bp, n - 1, system), system)), (double)n);
        if (sc_peek_error() != SC_OK && sc_peek_error() != saved) {
            sc_set_error(SC_ERR_CONVERGENCE);
            return b;
        }
        sc_complex d = sc_abs(sc_sub(nx, bp));
        b = nx;
        if (d.real < DBL_EPSILON && d.imag < DBL_EPSILON) break;
        bp = b;
    }
    return b;
}

sc_complex sc_root(sc_complex a, int n, double system) { return sc_root_iter(a, n, system, 100); }

static sc_complex sc_exp_iter(sc_complex a, double system, int max_iterations) {
    if (a.imag == 0.0) return sc_make(exp(a.real), 0.0);
    if (sc_norm_norm(a) > 1.0) {
        sc_complex h = sc_exp_iter(sc_shr(a, 1), system, max_iterations);
        return sc_mul(h, h, system);
    }
    sc_complex temp = sc_make(0.0, a.imag);
    sc_complex b = sc_one();
    for (int i = max_iterations; i >= 0; --i)
        b = sc_add_real(sc_div_real(sc_mul(temp, b, system), (double)(i + 1)), 1.0);
    return sc_mul_real(b, exp(a.real));
}

sc_complex sc_exp(sc_complex a, double system) { return sc_exp_iter(a, system, 20); }

sc_complex sc_polar(double rho, double theta, double system) {
    return sc_exp(sc_make(log(rho), theta), system);
}

static sc_complex sc_agm(sc_complex a, sc_complex b, double system, int max_iterations) {
    sc_complex s = sc_add(a, b);
    sc_complex p = sc_mul(a, b, system);
    for (int i = 0; i < max_iterations; ++i) {
        a = sc_div_real(s, 2.0);
        b = sc_sqrt(p, system);
        s = sc_add(a, b);
        p = sc_mul(a, b, system);
    }
    return a;
}

static sc_complex sc_log_iter(sc_complex a, double system, int max_iterations) {
    if (sc_is_zero(a)) { sc_set_error(SC_ERR_DOMAIN); return sc_zero(); }
    if (a.imag == 0.0 && a.real > 0.0) return sc_make(log(a.real), 0.0);
    if (sc_norm(a, system) < 0.0) { sc_set_error(SC_ERR_DOMAIN); return sc_zero(); }
    double L = sc_norm_norm(a);
    if (L > 1.0) return sc_add_real(sc_log_iter(sc_shr(a, 1), system, max_iterations), 0.69314718055994530942);
    if (L < 0.25) return sc_sub_real(sc_log_iter(sc_shl(a, 1), system, max_iterations), 0.69314718055994530942);
    sc_complex q = sc_one();
    for (int i = 1; i < max_iterations; ++i) q = sc_div_real(q, 2.0);
    sc_complex ans1 = sc_agm(sc_one(), q, system, 8);
    sc_complex ans2 = sc_agm(sc_one(), sc_mul(a, q, system), system, 8);
    if (sc_norm_norm(ans1) < DBL_EPSILON || sc_norm_norm(ans2) < DBL_EPSILON) {
        sc_set_error(SC_ERR_CONVERGENCE);
        return sc_zero();
    }
    if (sc_peek_error() != SC_OK) {
        sc_set_error(SC_ERR_CONVERGENCE);
        return sc_zero();
    }
    return sc_mul_real(sc_sub(sc_div(sc_one(), ans1, system), sc_div(sc_one(), ans2, system)), sc_pi_2());
}

sc_complex sc_log(sc_complex a, double system) { return sc_log_iter(a, system, 40); }

sc_complex sc_log_base(sc_complex a, sc_complex base, double system) {
    return sc_div(sc_log(a, system), sc_log(base, system), system);
}

double sc_arg(sc_complex a, double system) { return sc_log(a, system).imag; }

static sc_complex sc_cos_iter(sc_complex a, double system, int max_iterations) {
    if (a.imag == 0.0) return sc_make(cos(a.real), 0.0);
    while (a.real > sc_2pi()) a.real -= sc_2pi();
    while (a.real < -sc_2pi()) a.real += sc_2pi();
    if (sc_norm_norm(a) > 1.0) {
        sc_complex h = sc_cos_iter(sc_shr(a, 1), system, max_iterations);
        return sc_sub_real(sc_mul_real(sc_mul(h, h, system), 2.0), 1.0);
    }
    sc_complex b = sc_one();
    for (int i = max_iterations; i >= 0; --i) {
        b = sc_div_real(sc_mul(a, b, system), (double)(i + 1));
        if (i % 2) { /* odd: keep product only */ }
        else if (i % 4) b = sc_sub_real(b, 1.0);
        else b = sc_add_real(b, 1.0);
    }
    return b;
}

sc_complex sc_cos(sc_complex a, double system) { return sc_cos_iter(a, system, 20); }

sc_complex sc_sin(sc_complex a, double system) {
    return sc_cos(sc_sub(sc_make(sc_pi_2(), 0.0), a), system);
}

sc_complex sc_tan(sc_complex a, double system) { return sc_div(sc_sin(a, system), sc_cos(a, system), system); }

sc_complex sc_cot(sc_complex a, double system) { return sc_div(sc_cos(a, system), sc_sin(a, system), system); }

sc_complex sc_sec(sc_complex a, double system) { return sc_div(sc_one(), sc_cos(a, system), system); }

sc_complex sc_csc(sc_complex a, double system) { return sc_div(sc_one(), sc_sin(a, system), system); }

sc_complex sc_sinh(sc_complex a, double system) {
    return sc_div_real(sc_sub(sc_exp(a, system), sc_exp(sc_neg(a), system)), 2.0);
}

sc_complex sc_cosh(sc_complex a, double system) {
    return sc_div_real(sc_add(sc_exp(a, system), sc_exp(sc_neg(a), system)), 2.0);
}

sc_complex sc_tanh(sc_complex a, double system) { return sc_div(sc_sinh(a, system), sc_cosh(a, system), system); }

sc_complex sc_coth(sc_complex a, double system) { return sc_div(sc_cosh(a, system), sc_sinh(a, system), system); }

sc_complex sc_sech(sc_complex a, double system) { return sc_div(sc_one(), sc_cosh(a, system), system); }

sc_complex sc_csch(sc_complex a, double system) { return sc_div(sc_one(), sc_sinh(a, system), system); }

sc_complex sc_pow(sc_complex base, sc_complex exponent, double system) {
    if (sc_is_zero(exponent)) return sc_one();
    if (sc_is_zero(base)) return sc_zero();
    if (exponent.imag == 0.0 && exponent.real == round(exponent.real) &&
        exponent.real >= (double)INT_MIN && exponent.real <= (double)INT_MAX)
        return sc_pow_int(base, (int)exponent.real, system);
    return sc_exp(sc_mul(sc_log(base, system), exponent, system), system);
}

sc_complex sc_root_complex(sc_complex base, sc_complex exponent, double system) {
    if (sc_equal(exponent, sc_one())) return base;
    if (sc_is_zero(base)) return sc_zero();
    if (exponent.imag == 0.0 && exponent.real == round(exponent.real))
        return sc_root(base, (int)exponent.real, system);
    return sc_exp(sc_div(sc_log(base, system), exponent, system), system);
}

sc_complex sc_asinh(sc_complex a, double system) {
    return sc_log(sc_add(a, sc_sqrt(sc_add(sc_mul(a, a, system), sc_one()), system)), system);
}

sc_complex sc_acosh(sc_complex a, double system) {
    return sc_log(sc_add(a, sc_sqrt(sc_sub(sc_mul(a, a, system), sc_one()), system)), system);
}

sc_complex sc_atanh(sc_complex a, double system) {
    return sc_div_real(sc_log(sc_div(sc_add(sc_one(), a), sc_sub(sc_one(), a), system), system), 2.0);
}

sc_complex sc_acoth(sc_complex a, double system) {
    return sc_div_real(sc_log(sc_div(sc_add(a, sc_one()), sc_sub(a, sc_one()), system), system), 2.0);
}

sc_complex sc_asech(sc_complex a, double system) {
    return sc_log(sc_div(sc_add(sc_sqrt(sc_sub(sc_one(), sc_mul(a, a, system)), system), sc_one()), a, system), system);
}

sc_complex sc_acsch(sc_complex a, double system) {
    return sc_log(sc_div(sc_add(sc_sqrt(sc_add(sc_one(), sc_mul(a, a, system)), system), sc_one()), a, system), system);
}

sc_complex sc_integral(sc_complex a, sc_complex b, sc_func_s func, double system, double delta) {
    sc_complex d = sc_mul_real(sc_sub(b, a), delta);
    sc_complex ans = sc_zero();
    sc_complex ins = sc_mul(func(a, system), d, system);
    for (int i = 0; delta * (double)i <= 1.0; ++i) {
        ans = sc_add(ans, ins);
        a = sc_add(a, d);
        ins = sc_mul(func(a, system), d, system);
    }
    return ans;
}

sc_complex sc_integral_default(sc_complex a, sc_complex b, sc_func_s func, double system) {
    return sc_integral(a, b, func, system, 1e-6);
}

sc_complex sc_differential(sc_complex a, sc_func_s func, double system, double delta) {
    sc_complex d_real = sc_div_real(sc_sub(func(sc_add_real(a, delta), system), func(a, system)), delta);
    sc_complex di = sc_make(0.0, delta);
    sc_complex d_imag = sc_div(sc_sub(func(sc_add(a, di), system), func(a, system)), di, system);
    sc_complex d = sc_abs(sc_div(sc_sub(d_real, d_imag), func(a, system), system));
    if (d.real > 1e-3 || d.imag > 1e-3) sc_set_error(SC_ERR_DOMAIN);
    return sc_div_real(sc_add(d_real, d_imag), 2.0);
}

sc_complex sc_differential_default(sc_complex a, sc_func_s func, double system) {
    return sc_differential(a, func, system, 1e-6);
}

sc_complex sc_invert(sc_complex a, sc_func_s func, double system, int max_iterations) {
    sc_complex b = sc_one();
    sc_complex bp = b;
    for (int i = 0; i < max_iterations; ++i) {
        sc_error_t saved = sc_peek_error();
        sc_complex d1 = sc_differential(bp, func, system, 1e-6);
        sc_complex d2 = sc_differential(bp, func, system, 1e-6);
        if (sc_peek_error() != SC_OK && sc_peek_error() != saved) {
            sc_set_error(SC_ERR_CONVERGENCE);
            return b;
        }
        b = sc_add(sc_sub(bp, sc_div(func(bp, system), d1, system)), sc_div(a, d2, system));
        sc_complex d = sc_abs(sc_sub(b, bp));
        if (d.real < 1e-3 && d.imag < 1e-3) break;
        bp = b;
    }
    return b;
}

sc_complex sc_invert_default(sc_complex a, sc_func_s func, double system) {
    return sc_invert(a, func, system, 100);
}

sc_complex sc_acos(sc_complex a, double system) { return sc_invert(a, sc_cos, system, 100); }

sc_complex sc_asin(sc_complex a, double system) {
    return sc_sub(sc_make(sc_pi_2(), 0.0), sc_acos(a, system));
}

sc_complex sc_atan(sc_complex a, double system) {
    return sc_acos(sc_div(sc_one(), sc_sqrt(sc_add(sc_pow_int(a, 2, system), sc_one()), system), system), system);
}

sc_complex sc_acot(sc_complex a, double system) {
    return sc_sub(sc_make(sc_pi_2(), 0.0), sc_atan(a, system));
}

sc_complex sc_asec(sc_complex a, double system) { return sc_acos(sc_div(sc_one(), a, system), system); }

sc_complex sc_acsc(sc_complex a, double system) {
    return sc_sub(sc_make(sc_pi_2(), 0.0), sc_asec(a, system));
}

static sc_complex sc_gamma_iter(sc_complex a, double system, int max_iterations) {
    if (sc_equal(a, sc_one())) return sc_one();
    if (a.imag == 0.0 && a.real < 0.0 && round(a.real) == a.real) {
        sc_set_error(SC_ERR_DOMAIN);
        return sc_zero();
    }
    if (a.real < 0.5) {
        return sc_real_div(sc_pi(), sc_mul(sc_sin(sc_mul_real(a, sc_pi()), system),
                                           sc_gamma_iter(sc_sub(sc_one(), a), system, max_iterations), system), system);
    }
    if (a.real > 1.0)
        return sc_mul(sc_gamma_iter(sc_sub_real(a, 1.0), system, max_iterations), sc_sub_real(a, 1.0), system);
    sc_complex b = sc_neg(sc_log(a, system));
    for (int i = 1; i <= max_iterations; ++i) {
        b = sc_sub(b, sc_log(sc_add_real(sc_div_real(a, (double)i), 1.0), system));
        b = sc_add(b, sc_mul_real(a, log(1.0 + 1.0 / i)));
    }
    return sc_exp(b, system);
}

sc_complex sc_gamma(sc_complex a, double system) { return sc_gamma_iter(a, system, 100000); }

sc_complex sc_beta(sc_complex a, sc_complex b, double system) {
    return sc_div(sc_mul(sc_gamma(a, system), sc_gamma(b, system), system),
                  sc_gamma(sc_add(a, b), system), system);
}

sc_complex sc_square(sc_complex x, double system) { return sc_mul(x, x, system); }

sc_complex sc_M(sc_complex a, sc_complex c, sc_func_s func, double system, bool L) {
    sc_complex v = L ? sc_sub(a, func(a, system)) : func(a, system);
    return sc_add(v, c);
}

sc_complex sc_T(sc_complex a, sc_complex c, sc_func_s func, double system, bool L) {
    return sc_M(sc_conj(a), c, func, system, L);
}

sc_complex sc_B(sc_complex a, sc_complex c, sc_func_s func, double system, bool L) {
    return sc_M(sc_abs(a), c, func, system, L);
}

sc_complex sc_hyper_xexp(sc_complex a, int n, double system) {
    sc_complex ans = a;
    for (int i = 1; i <= n; ++i) ans = sc_exp(ans, system);
    for (int i = 1; i <= -n; ++i) ans = sc_log(ans, system);
    return sc_mul(ans, a, system);
}

static int sc_hyper_internal_n = 0;

static double sc_hyper_internal_system = -1.0;

static sc_complex sc_hyper_xexp_n(sc_complex x, double system) {
    (void)system;
    return sc_hyper_xexp(x, sc_hyper_internal_n, sc_hyper_internal_system);
}

sc_complex sc_hyper_omega(sc_complex a, int n, double system) {
    if (n >= 0) {
        sc_hyper_internal_n = n;
        sc_hyper_internal_system = system;
        return sc_invert(a, sc_hyper_xexp_n, system, 100);
    }
    sc_complex ans = sc_hyper_omega(a, -n, system);
    for (int i = 1; i <= -n; ++i) ans = sc_exp(ans, system);
    return ans;
}
