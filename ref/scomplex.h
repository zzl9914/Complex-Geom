/**
 * @file scomplex.h
 * @brief Complex / split-complex / dual-number arithmetic (C API).
 *
 * Semantics match this tree's C++ port in scomplex.hpp (and the original
 * scomplex_no_throw.hpp algebra). Errors are reported, not swallowed:
 * failing paths write sc_last_error.
 *
 * system = -1 ordinary complex (i^2 = -1),
 *          +1 split / hyperbolic (j^2 = +1),
 *           0 dual (epsilon^2 = 0).
 *
 * Public symbols use the sc_ prefix. Implementations live in scomplex.c;
 * compile and link that translation unit with the caller.
 */

#ifndef SCOMPLEX_H_
#define SCOMPLEX_H_

#include <stdbool.h>
#include <stdint.h>

#ifdef __cplusplus
extern "C" {
#endif

#define SC_SYSTEM_COMPLEX (-1.0)
#define SC_SYSTEM_SPLIT   ( 1.0)
#define SC_SYSTEM_DUAL    ( 0.0)

typedef struct {
    double real;
    double imag;
} sc_complex;

typedef sc_complex (*sc_func_s)(sc_complex, double system);

typedef enum {
    SC_OK = 0,
    SC_ERR_DIV_BY_ZERO,
    SC_ERR_OVERFLOW,
    SC_ERR_DOMAIN,
    SC_ERR_CONVERGENCE,
    SC_ERR_INVALID
} sc_error_t;

extern sc_error_t sc_last_error;

const char* sc_error_string(sc_error_t e);
void sc_set_error(sc_error_t e);
sc_error_t sc_get_error(void);
sc_error_t sc_peek_error(void);
float sc_flnot(float a);
float sc_flor(float a, float b);
float sc_fland(float a, float b);
float sc_flxor(float a, float b);
double sc_dblnot(double a);
double sc_dblor(double a, double b);
double sc_dbland(double a, double b);
double sc_dblxor(double a, double b);
double sc_shift_left(double a, uint64_t b);
double sc_shift_right(double a, uint64_t b);
double sc_round_left(double a);
double sc_pi(void);
double sc_2pi(void);
double sc_pi_2(void);
double sc_pi_3(void);
double sc_e(void);
double sc_sqrt2(void);
double sc_sqrt1_2(void);
double sc_cbrt2(void);
double sc_cbrt1_2(void);
sc_complex sc_make(double real, double imag);
sc_complex sc_make_real(double real);
sc_complex sc_zero(void);
sc_complex sc_one(void);
sc_complex sc_i(void);
sc_complex sc_exp_i_pi_1_6(void);
sc_complex sc_exp_i_pi_1_9(void);
sc_complex sc_from_string(const char* str);
double sc_real(sc_complex a);
double sc_imag(sc_complex a);
bool sc_equal(sc_complex a, sc_complex b);
bool sc_not_equal(sc_complex a, sc_complex b);
bool sc_less(sc_complex a, sc_complex b);
bool sc_greater(sc_complex a, sc_complex b);
bool sc_less_equal(sc_complex a, sc_complex b);
bool sc_greater_equal(sc_complex a, sc_complex b);
bool sc_is_zero(sc_complex a);
bool sc_nonzero(sc_complex a);
bool sc_is_finite(sc_complex a);
bool sc_is_nan(sc_complex a);
sc_complex sc_add(sc_complex a, sc_complex b);
sc_complex sc_sub(sc_complex a, sc_complex b);
sc_complex sc_mul(sc_complex a, sc_complex b, double system);
sc_complex sc_div(sc_complex a, sc_complex b, double system);
sc_complex sc_mod(sc_complex a, sc_complex b, double system);
sc_complex sc_neg(sc_complex a);
sc_complex sc_conj(sc_complex a);
sc_complex sc_add_real(sc_complex a, double b);
sc_complex sc_sub_real(sc_complex a, double b);
sc_complex sc_real_sub(double a, sc_complex b);
sc_complex sc_mul_real(sc_complex a, double b);
sc_complex sc_div_real(sc_complex a, double b);
sc_complex sc_real_div(double a, sc_complex b, double system);
sc_complex sc_mod_real(sc_complex a, double b);
sc_complex sc_real_mod(double a, sc_complex b, double system);
sc_complex sc_shl(sc_complex a, uint64_t b);
sc_complex sc_shr(sc_complex a, uint64_t b);
double sc_norm(sc_complex a, double system);
double sc_norm_norm(sc_complex a);
sc_complex sc_abs(sc_complex a);
sc_complex sc_fabs(sc_complex a);
double sc_abs2(sc_complex a);
double sc_modulus(sc_complex a);
sc_complex sc_pow_int(sc_complex base, int exponent, double system);
sc_complex sc_sqrt(sc_complex a, double system);
sc_complex sc_cbrt(sc_complex a, double system);
sc_complex sc_root(sc_complex a, int n, double system);
sc_complex sc_exp(sc_complex a, double system);
sc_complex sc_polar(double rho, double theta, double system);
sc_complex sc_log(sc_complex a, double system);
sc_complex sc_log_base(sc_complex a, sc_complex base, double system);
double sc_arg(sc_complex a, double system);
sc_complex sc_cos(sc_complex a, double system);
sc_complex sc_sin(sc_complex a, double system);
sc_complex sc_tan(sc_complex a, double system);
sc_complex sc_cot(sc_complex a, double system);
sc_complex sc_sec(sc_complex a, double system);
sc_complex sc_csc(sc_complex a, double system);
sc_complex sc_sinh(sc_complex a, double system);
sc_complex sc_cosh(sc_complex a, double system);
sc_complex sc_tanh(sc_complex a, double system);
sc_complex sc_coth(sc_complex a, double system);
sc_complex sc_sech(sc_complex a, double system);
sc_complex sc_csch(sc_complex a, double system);
sc_complex sc_pow(sc_complex base, sc_complex exponent, double system);
sc_complex sc_root_complex(sc_complex base, sc_complex exponent, double system);
sc_complex sc_asinh(sc_complex a, double system);
sc_complex sc_acosh(sc_complex a, double system);
sc_complex sc_atanh(sc_complex a, double system);
sc_complex sc_acoth(sc_complex a, double system);
sc_complex sc_asech(sc_complex a, double system);
sc_complex sc_acsch(sc_complex a, double system);
sc_complex sc_integral(sc_complex a, sc_complex b, sc_func_s func, double system, double delta);
sc_complex sc_integral_default(sc_complex a, sc_complex b, sc_func_s func, double system);
sc_complex sc_differential(sc_complex a, sc_func_s func, double system, double delta);
sc_complex sc_differential_default(sc_complex a, sc_func_s func, double system);
sc_complex sc_invert(sc_complex a, sc_func_s func, double system, int max_iterations);
sc_complex sc_invert_default(sc_complex a, sc_func_s func, double system);
sc_complex sc_acos(sc_complex a, double system);
sc_complex sc_asin(sc_complex a, double system);
sc_complex sc_atan(sc_complex a, double system);
sc_complex sc_acot(sc_complex a, double system);
sc_complex sc_asec(sc_complex a, double system);
sc_complex sc_acsc(sc_complex a, double system);
sc_complex sc_gamma(sc_complex a, double system);
sc_complex sc_beta(sc_complex a, sc_complex b, double system);
sc_complex sc_square(sc_complex x, double system);
sc_complex sc_M(sc_complex a, sc_complex c, sc_func_s func, double system, bool L);
sc_complex sc_T(sc_complex a, sc_complex c, sc_func_s func, double system, bool L);
sc_complex sc_B(sc_complex a, sc_complex c, sc_func_s func, double system, bool L);
sc_complex sc_hyper_xexp(sc_complex a, int n, double system);
sc_complex sc_hyper_omega(sc_complex a, int n, double system);

#ifdef __cplusplus
}
#endif

#endif /* SCOMPLEX_H_ */
