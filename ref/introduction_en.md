# scomplex

Other languages: [使用说明_中文.md](使用说明_中文.md) · [Mode_d'emploi_français.md](Mode_d'emploi_français.md) · [使用説明_日本語.md](使用説明_日本語.md)

scomplex is a small library for arithmetic in three related two-dimensional
algebras: ordinary complex numbers, split-complex (hyperbolic) numbers, and
dual numbers. The same multiplication law is used everywhere; a scalar
parameter `system` chooses which algebra you are in.

There are three language ports of one algebra:

| Language | Files | How the algebra is selected |
|---|---|---|
| C | `scomplex.h`, `scomplex.c` | explicit `double system` argument |
| C++ | `scomplex.hpp`, `scomplex.cpp` | `complex<_Tp>::set_system` / `get_system` |
| Haskell | `Scomplex.hs` (module `SComplex`) | explicit `Double` argument, first parameter of most functions |

The ports are meant to agree with each other and with the original C++ source
`scomplex_no_throw.hpp`. That original lived as a “written in the standard
library” overlay: it sits in `std::complex::experimental::zhangzl`, piles on
ordinary standard headers, and **does not** `#include <complex>`, so it does
not collide with `std::complex`. The C++ port here keeps that overlay. The
one intentional departure from the canvas `no_throw` build is error handling:
failing paths actually report errors instead of catching and swallowing them.

Transcendental functions are **not** `std::sin` / `atan2` / `std::polar`.
Those formulas are valid only for ordinary complex numbers (`i² = −1`).
Here they are implemented with the library’s own multiplication (power
series, Newton iteration, arithmetic–geometric mean).

---

## Files

```
scomplex.h            C declarations, types, sc_last_error
scomplex.c            C implementations (internal helpers are file-static)
scomplex.hpp          C++ class template, free function templates, exceptions
scomplex.cpp          non-template C++ helpers that the header calls
Scomplex.hs           Haskell module SComplex
introduction_en.md           this file (English)
使用说明_中文.md             Chinese
Mode_d'emploi_français.md    French
使用説明_日本語.md           Japanese
```

`scomplex.cpp.h` was a misnamed header and is gone. C++ callers include
`scomplex.hpp` and **link** `scomplex.cpp`.

Templates that take a user functor (`differential`, `invert`, `integral`)
cannot move into `scomplex.cpp`: the functor type is only known in the
caller’s translation unit. The class template and the other free templates
stay in the header for the same reason.

---

## The three algebras

A value is always a pair `(x, y)` of real components. Multiplication is

```
(x, y) * (u, v) = (x u + system * y v,  x v + y u)
```

| `system` | Algebra | Imaginary unit | Square |
|---|---|---|---|
| `−1` | ordinary complex | `i` | `i² = −1` |
| `+1` | split / hyperbolic | `j` | `j² = +1` |
| `0` | dual | `ε` | `ε² = 0` |

C macros: `SC_SYSTEM_COMPLEX`, `SC_SYSTEM_SPLIT`, `SC_SYSTEM_DUAL`.

C++ default is `complex<double>::_system = −1`. Change it with
`complex<double>::set_system(+1)` or `set_system(0)` before operations that
read the static. Haskell’s `defaultSystem` is `−1`.

The (algebraic) norm used for domain checks is

```
norm(z) = x² − system * y²
```

`norm_norm(z) = x² + y²` is the Euclidean square length. It does **not**
depend on `system`. Scaling by powers of two (`shift_left` / `shift_right`,
C `sc_shl` / `sc_shr`) uses `norm_norm` to keep Newton / AGM arguments in a
comfortable range.

---

## Error handling

Errors are real. They are not discarded.

| Port | Mechanism |
|---|---|
| C | `sc_last_error` (`SC_OK`, `SC_ERR_DIV_BY_ZERO`, `SC_ERR_OVERFLOW`, `SC_ERR_DOMAIN`, `SC_ERR_CONVERGENCE`, `SC_ERR_INVALID`). `sc_set_error` writes it. `sc_peek_error` reads it. `sc_get_error` reads and clears. `sc_error_string` maps a code to text. |
| C++ | `throw`. Custom types `complex_error`, `complex_division_by_zero`, `complex_overflow`, `complex_domain_error`, `complex_convergence_error` inherit `std::runtime_error`. Some original paths still throw `std::domain_error` or `std::invalid_argument` (string constructor, `sqrt` domain). |
| Haskell | `throw` of `ComplexError` (`DivisionByZero`, `Overflow`, `DomainError String`, `ConvergenceError`, `InvalidFormat String`). |

`sc_last_error` is one process-wide object (defined in `scomplex.c`). C++
`_system` and `print_number` are static members of the class template, so
they are per-type, process-wide, and not thread-safe. Do not call the C or
C++ libraries concurrently on the same type without external locking.

Overflow is checked when a pair is constructed: non-finite components, or
magnitude above `max/2` (C: `DBL_MAX / 2`). Division by zero is
`|denominator| < ε`.

---

## Arithmetic and comparisons

Addition, subtraction, negation, and conjugation do not use `system`.
Multiplication, division, and remainder do.

C names: `sc_add`, `sc_sub`, `sc_mul`, `sc_div`, `sc_mod`, `sc_neg`,
`sc_conj`, plus mixed real helpers (`sc_add_real`, `sc_mul_real`,
`sc_div_real`, `sc_real_div`, …).

C++ uses ordinary operators on `complex<_Tp>`, including `+=`, `*=` and so
on. Remainder uses `round_left` on the quotient (fractional part toward
zero-ish via `x - round(x)`), then multiplies back.

Comparisons are lexicographic: real first, then imag. C: `sc_equal`,
`sc_less`, … Haskell: derived `Eq` only; there is no `Ord` instance.

**`abs` is component-wise**, `|x| + |y| i`, not the modulus. Euclidean
length is `sc_modulus` / `sqrt(norm_norm(z))`. This matches the original.

---

## Bit operations and exponent shifts

IEEE-754 bit tricks exist because the original used them as cheap
power-of-two scaling:

- C: `sc_flnot` / `sc_flor` / `sc_fland` / `sc_flxor` (float),
  `sc_dblnot` / `sc_dblor` / `sc_dbland` / `sc_dblxor` (double),
  `sc_shift_left` / `sc_shift_right` (add/subtract from the exponent field).
- C++: `digital_not` / `digital_or` / `digital_and` / `digital_xor`,
  `shift_left` / `shift_right`, `round_left`, `medium_mod`.
- Haskell: `shiftLeft` / `shiftRight` via `scaleFloat`.

`shift_left(a, n)` is `a * 2ⁿ` on a finite non-zero double, implemented by
adding `n` to the exponent bits. Zero stays zero. These are **not** C
integer shifts.

---

## Parsing

Accepted forms (whitespace stripped):

- real only: `3`, `-2.5`
- imaginary only: `i`, `+i`, `-i`
- both: `3+4i`, `3-4i`, `3+i`

C++ uses a regex in the `std::string` constructor and throws
`std::invalid_argument` on mismatch. C `sc_from_string` writes
`SC_ERR_INVALID` and returns zero. Haskell `fromString` throws
`InvalidFormat`.

---

## Elementary functions

All of these take `system` (C/Haskell) or read `_system` (C++).

### Square root, cube root, integer root

Newton iteration after reducing `|z|` into `[0.5, 2]` by exponent shifts.

- Real non-negative square roots (and real cube roots) use the host `sqrt` /
  `cbrt`.
- If `norm(z) < 0`, **square root is a domain error** (ordinary negative
  reals have `norm = x² ≥ 0`, so they are allowed; split `sqrt(j)` is not).
- Dual `sqrt(−1)` has `norm = 1`, so it is **not** a domain error. The
  iteration can produce a huge imaginary part. That is original behaviour.

Initial guesses: `1` in the right half-plane, `exp(i π / 6)` (or conjugate)
otherwise for square root; `exp(i π / 9)` for cube root.

### Exponential

If the argument is real, host `exp`. If `norm_norm(z) > 1`, use
`exp(z) = exp(z/2)²` (with library multiply). Otherwise a Horner series in
the imaginary part, then multiply by `exp(real)`. Default 20 terms.

Consequence: `exp(i π) ≈ −1` in the ordinary system, `exp(j) = cosh 1 + j sinh 1`
in the split system, `exp(ε) = 1 + ε` in the dual system.

`polar(ρ, θ)` is `exp(log ρ + θ · unit)`, using this `exp`, not `std::polar`.

### Logarithm and argument

If `z = 0` or `norm(z) < 0`, domain error. Positive reals use host `log`.
Otherwise scale into a band with exponent shifts (adding `ln 2`), then an
AGM identity:

```
log(z) ≈ (π/2) * (1/AGM(1, q) − 1/AGM(1, z q))
```

with a tiny `q`. Default 40 outer steps, 8 AGM steps. `arg(z)` is
`imag(log(z))`.

### Cosine and sine

Cosine: reduce the real part modulo `2π`, halve until `norm_norm ≤ 1`
(`cos 2α = 2 cos²α − 1`), then a series. Default 20 terms.

**Sine is `cos(π/2 − z)`**, as in the original, not an independent series.

`tan`, `cot`, `sec`, `csc` are ratios of those.

### Hyperbolic functions

```
sinh z = (exp(z) − exp(−z)) / 2
cosh z = (exp(z) + exp(−z)) / 2
```

and the usual ratios. Inverse hyperbolics are the usual log/sqrt identities,
evaluated with this library’s `log` and `sqrt`.

### Power

Integer real exponents use binary exponentiation (`sc_pow_int`). Otherwise
`exp(log(base) * exponent)`. Complex roots are `exp(log(base) / exponent)`,
or Newton `root` when the exponent is an integer.

---

## Inverse trigonometric functions

There are no closed-form `asin` series. `acos` is Newton inversion of `cos`
(`invert` / `sc_invert`). The others are identities on that:

```
asin(z) = π/2 − acos(z)
atan(z) = acos(1 / sqrt(z² + 1))
acot(z) = π/2 − atan(z)
asec(z) = acos(1/z)
acsc(z) = π/2 − asec(z)
```

**`atan(−1) = +π/4`**, not `−π/4`. The identity drops the sign on the
negative real axis. That is original, not a port bug.

`invert` differentiates the target function twice per step (two finite
differences) and stops when the component-wise `abs` of the update is below
`10⁻³`, or after `max_iterations` (default 100). If a differentiation
throws, C++ rethrows `complex_convergence_error`.

---

## Gamma and beta

Weierstrass product, default **100 000** iterations (original). Reflection
for `Re(z) < 1/2`, recurrence for `Re(z) > 1`. Non-positive integers on the
real axis are a domain error. `beta(a, b) = Γ(a) Γ(b) / Γ(a+b)`.

---

## Calculus

Finite-difference derivative at `a`:

```
d_real = (f(a + δ) − f(a)) / δ
d_imag = (f(a + i δ) − f(a)) / (i δ)
```

If those two disagree (component-wise `abs` of the relative gap `> 10⁻³`),
the function is treated as not complex-differentiable: C writes
`SC_ERR_DOMAIN`, C++/Haskell throw. Otherwise the result is the average.
Default `δ = 10⁻⁶`.

Integral from `a` to `b` is a left Riemann sum with step `δ (b − a)`,
running while `δ · k ≤ 1`. Default `δ = 10⁻⁶`.

C callbacks have type `sc_func_s`: `sc_complex (*)(sc_complex, double system)`.
C++ and Haskell take any compatible functor / function.

---

## Hyperoperations and M / T / B

`hyper_xexp(z, n)` applies `exp` `n` times (or `log` `|n|` times if `n < 0`)
and then multiplies by `z`.

`hyper_omega(z, n)` is the inverse of that in `z` (Newton `invert`) for
`n ≥ 0`; negative `n` applies extra `exp`s to `hyper_omega(z, −n)`.

C keeps the current `n` and `system` in file-static variables while
`sc_hyper_omega` calls `sc_invert`. That is another reason the C API is not
re-entrant.

`M`, `T`, `B` (and C++ batch overloads on arrays) are the original
Mandelbrot-style maps: `M(a, c, f, L)` is `a − f(a) + c` or `f(a) + c`
according to `L`; `T` conjugates first; `B` takes component-wise `abs`
first. `sc_square` / C++ default `pow(x, 2)` is the usual quadratic.

---

## C++ overlay details

```cpp
#include "scomplex.hpp"
using namespace std::complex::experimental::zhangzl;

complex<double>::set_system(-1.0);
complex<double> z(2.0, 1.0);
auto w = exp(z);
complex<double> i("i");
```

- Namespace: `std::complex::experimental::zhangzl`.
- Do **not** `#include <complex>` in the same way you would for
  `std::complex`; this header already defines `complex` in a nested
  `std::complex` namespace.
- Compile with `/utf-8` on MSVC. The sources are UTF-8; without that flag
  MSVC may emit C4819 and can eat the definition of `_system` (LNK2019).
- Link `scomplex.cpp`. Without it, `shift_left` / `shift_right` /
  `round_left` / `medium_mod` are unresolved.
- `COMPLEX_ALWAYS_INLINE` is `__forceinline` on MSVC. Combined with
  `inline` on some operators this yields C4141; it does not affect linking.
- `sqrt` of a negative-norm value throws `std::domain_error` with
  `"Square root of non-positive number"`, as in the original, not
  `complex_domain_error`.

---

## C details

```c
#include "scomplex.h"

sc_complex a = sc_make(2.0, 1.0);
sc_complex b = sc_make(1.0, 3.0);
sc_complex p = sc_mul(a, b, SC_SYSTEM_COMPLEX); /* -1 + 7i */
sc_error_t e = sc_get_error();
```

Compile and link `scomplex.c`. C11 is enough. `extern "C"` is provided for
C++ callers of the C API.

`sc_get_error` clears the flag. If you need to inspect without consuming,
use `sc_peek_error`.

---

## Haskell details

```haskell
import qualified SComplex as SC

let z = SC.fromRealImag 2 1
    w = SC.multiply SC.defaultSystem z (SC.fromRealImag 1 3)
```

Most operations take `system` first. `add` / `subtract` / `conjugate` /
`fromString` / `abs` do not, because they are system-independent.

Prelude `sqrt`, `exp`, `log`, `sin`, … are hidden; use `SComplex` or
`qualified Prelude as P`.

---

## Constants

C functions `sc_pi`, `sc_2pi`, `sc_pi_2`, `sc_pi_3`, `sc_e`, `sc_sqrt2`,
`sc_sqrt1_2`, `sc_cbrt2`, `sc_cbrt1_2`. C++ statics on `complex<_Tp>`
(`pi()`, `e()`, …) plus `i()`. Haskell locals `piVal`, `sqrt2`, …

`sc_exp_i_pi_1_6` / `sc_exp_i_pi_1_9` are the Newton seeds
`exp(i π / 6)` and `exp(i π / 9)`.

---

## Original semantics to keep

These look like textbook bugs if you score against principal values. They
are the original library:

1. `abs` is per-component, not modulus.
2. `sin(z) = cos(π/2 − z)`.
3. `atan(−1) = +π/4`.
4. Dual `sqrt(−1)` is not a domain error (`norm = 1`); the value can be huge.
5. Split `sqrt(j)` **is** a domain error (`norm = −1`).
6. Gamma uses 100 000 Weierstrass terms.
7. Finite-difference invert uses a `10⁻³` stop, so inverse trig is coarse.

---

## Building (MSVC / GHC)

Sources are UTF-8.

```text
cl /utf-8 /std:c11 scomplex.c your.c
cl /utf-8 /std:c++17 /EHsc scomplex.cpp your.cpp
ghc -i. YourMain.hs
```

GCC/Clang equivalents: compile `scomplex.c` or `scomplex.cpp` together with
the caller; C++ needs exceptions (`-fexceptions`, usually on by default).
The C++ header uses MSVC/GCC alignment and prefetch macros; it does not
require SIMD at run time.

---

## Suggested reading order in the sources

1. Multiplication and `norm` — the whole design is that one product.
2. `exp` / `log` / `sqrt` — how `system` actually changes analysis.
3. `cos` then `sin` — the reduction identity.
4. `invert` then `acos` / `atan` — why inverse trig looks the way it does.
5. Error paths in `scomplex.c` vs throws in `scomplex.hpp` vs `throwError`
   in `Scomplex.hs`.
