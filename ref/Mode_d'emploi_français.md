# scomplex

Autres langues : [introduction_en.md](introduction_en.md) · [使用说明_中文.md](使用说明_中文.md) · [使用説明_日本語.md](使用説明_日本語.md)

scomplex est une petite bibliothèque d’arithmétique sur trois algèbres
bidimensionnelles voisines : les nombres complexes ordinaires, les nombres
split-complexes (hyperboliques) et les nombres duaux. La même loi de
multiplication s’applique partout ; un scalaire `system` choisit l’algèbre.

Il existe trois ports d’une seule algèbre :

| Langage | Fichiers | Choix de l’algèbre |
|---|---|---|
| C | `scomplex.h`, `scomplex.c` | argument explicite `double system` |
| C++ | `scomplex.hpp`, `scomplex.cpp` | `complex<_Tp>::set_system` / `get_system` |
| Haskell | `Scomplex.hs` (module `SComplex`) | argument `Double` explicite, premier paramètre de la plupart des fonctions |

Les ports doivent s’accorder entre eux et avec la source C++ d’origine
`scomplex_no_throw.hpp`. Celle-ci se présentait comme une surcouche « écrite
dans la bibliothèque standard » : elle vit dans
`std::complex::experimental::zhangzl`, empile les en-têtes standard
ordinaires, et **n’** `#include` **pas** `<complex>`, afin de ne pas
entrer en collision avec `std::complex`. Le port C++ conserve cette
surcouche. L’écart volontaire par rapport à la construction `no_throw` du
canevas est la gestion d’erreurs : les chemins d’échec signalent vraiment
l’erreur au lieu de l’avaler après un `catch`.

Les fonctions transcendantes **ne sont pas** `std::sin` / `atan2` /
`std::polar`. Ces formules ne valent que pour les complexes ordinaires
(`i² = −1`). Ici elles sont réalisées avec la multiplication de la
bibliothèque (séries entières, Newton, moyenne arithmético-géométrique).

---

## Fichiers

```
scomplex.h            déclarations C, types, sc_last_error
scomplex.c            implémentations C (aides internes static au fichier)
scomplex.hpp          patron de classe C++, patrons de fonctions libres, exceptions
scomplex.cpp          aides C++ non-patron appelées par l’en-tête
Scomplex.hs           module Haskell SComplex
introduction_en.md           anglais
使用说明_中文.md             chinois
Mode_d'emploi_français.md    ce fichier (français)
使用説明_日本語.md           japonais
```

`scomplex.cpp.h` était un en-tête mal nommé ; il a disparu. Les appelants
C++ incluent `scomplex.hpp` et **lient** `scomplex.cpp`.

Les patrons qui prennent un foncteur utilisateur (`differential`,
`invert`, `integral`) ne peuvent pas migrer dans `scomplex.cpp` : le type
du foncteur n’est connu que dans l’unité de traduction de l’appelant. Le
patron de classe et les autres patrons de fonctions libres restent dans
l’en-tête pour la même raison.

---

## Les trois algèbres

Une valeur est toujours un couple `(x, y)` de composantes réelles. La
multiplication est

```
(x, y) * (u, v) = (x u + system * y v,  x v + y u)
```

| `system` | Algèbre | Unité imaginaire | Carré |
|---|---|---|---|
| `−1` | complexe ordinaire | `i` | `i² = −1` |
| `+1` | split / hyperbolique | `j` | `j² = +1` |
| `0` | dual | `ε` | `ε² = 0` |

Macros C : `SC_SYSTEM_COMPLEX`, `SC_SYSTEM_SPLIT`, `SC_SYSTEM_DUAL`.

Par défaut en C++ : `complex<double>::_system = −1`. Changer avec
`complex<double>::set_system(+1)` ou `set_system(0)` avant les opérations
qui lisent la statique. En Haskell, `defaultSystem` vaut `−1`.

La norme (algébrique) utilisée pour les contrôles de domaine est

```
norm(z) = x² − system * y²
```

`norm_norm(z) = x² + y²` est le carré de la longueur euclidienne. Elle
**ne** dépend **pas** de `system`. Les mises à l’échelle par puissances de
deux (`shift_left` / `shift_right`, `sc_shl` / `sc_shr` en C) s’appuient
sur `norm_norm` pour garder les arguments de Newton / AGM dans une plage
confortable.

---

## Gestion des erreurs

Les erreurs sont réelles. Elles ne sont pas jetées sous le tapis.

| Port | Mécanisme |
|---|---|
| C | `sc_last_error` (`SC_OK`, `SC_ERR_DIV_BY_ZERO`, `SC_ERR_OVERFLOW`, `SC_ERR_DOMAIN`, `SC_ERR_CONVERGENCE`, `SC_ERR_INVALID`). `sc_set_error` écrit. `sc_peek_error` lit. `sc_get_error` lit et efface. `sc_error_string` convertit un code en texte. |
| C++ | `throw`. Types personnalisés `complex_error`, `complex_division_by_zero`, `complex_overflow`, `complex_domain_error`, `complex_convergence_error`, dérivés de `std::runtime_error`. Certains chemins d’origine lancent encore `std::domain_error` ou `std::invalid_argument` (constructeur chaîne, domaine de `sqrt`). |
| Haskell | `throw` de `ComplexError` (`DivisionByZero`, `Overflow`, `DomainError String`, `ConvergenceError`, `InvalidFormat String`). |

`sc_last_error` est un objet unique au processus (défini dans `scomplex.c`).
`_system` et `print_number` en C++ sont des membres statiques du patron de
classe, donc une copie par type, globale au processus, et non sûrs vis-à-vis
des fils d’exécution. Ne pas appeler les bibliothèques C ou C++ en parallèle
sur le même type sans verrou externe.

Le dépassement est contrôlé à la construction d’un couple : composantes non
finies, ou magnitude au-delà de `max/2` (C : `DBL_MAX / 2`). La division par
zéro est `|dénominateur| < ε`.

---

## Arithmétique et comparaisons

Addition, soustraction, négation et conjugaison n’utilisent pas `system`.
Multiplication, division et reste, si.

Noms C : `sc_add`, `sc_sub`, `sc_mul`, `sc_div`, `sc_mod`, `sc_neg`,
`sc_conj`, plus les aides mixtes avec un réel (`sc_add_real`,
`sc_mul_real`, `sc_div_real`, `sc_real_div`, …).

Le C++ utilise les opérateurs ordinaires sur `complex<_Tp>`, y compris
`+=`, `*=`, etc. Le reste applique `round_left` au quotient (partie
fractionnaire via `x - round(x)`), puis remultiplie.

Les comparaisons sont lexicographiques : réel d’abord, puis imaginaire. C :
`sc_equal`, `sc_less`, … Haskell : seulement `Eq` dérivé ; pas d’instance
`Ord`.

**`abs` est composante par composante**, `|x| + |y| i`, pas le module. La
longueur euclidienne est `sc_modulus` / `sqrt(norm_norm(z))`. C’est
conforme à l’original.

---

## Opérations bit à bit et décalages d’exposant

Les tours IEEE-754 existent parce que l’original s’en servait comme mise à
l’échelle bon marché par puissances de deux :

- C : `sc_flnot` / `sc_flor` / `sc_fland` / `sc_flxor` (float),
  `sc_dblnot` / `sc_dblor` / `sc_dbland` / `sc_dblxor` (double),
  `sc_shift_left` / `sc_shift_right` (addition/soustraction dans le champ
  exposant).
- C++ : `digital_not` / `digital_or` / `digital_and` / `digital_xor`,
  `shift_left` / `shift_right`, `round_left`, `medium_mod`.
- Haskell : `shiftLeft` / `shiftRight` via `scaleFloat`.

`shift_left(a, n)` vaut `a * 2ⁿ` pour un double fini non nul, en ajoutant
`n` aux bits d’exposant. Zéro reste zéro. Ce **ne sont pas** des décalages
entiers du C.

---

## Analyse syntaxique

Formes acceptées (blancs retirés) :

- réel seul : `3`, `-2.5`
- imaginaire seul : `i`, `+i`, `-i`
- les deux : `3+4i`, `3-4i`, `3+i`

Le C++ utilise une expression rationnelle dans le constructeur
`std::string` et lance `std::invalid_argument` en cas d’échec. En C,
`sc_from_string` écrit `SC_ERR_INVALID` et renvoie zéro. En Haskell,
`fromString` lance `InvalidFormat`.

---

## Fonctions élémentaires

Toutes prennent `system` (C/Haskell) ou lisent `_system` (C++).

### Racine carrée, cubique, entière

Itération de Newton après réduction de `|z|` dans `[0.5, 2]` par décalages
d’exposant.

- Les racines carrées réelles non négatives (et les racines cubiques réelles)
  utilisent `sqrt` / `cbrt` de l’hôte.
- Si `norm(z) < 0`, **la racine carrée est une erreur de domaine** (un réel
  négatif ordinaire a `norm = x² ≥ 0`, donc il est permis ; `sqrt(j)` split
  ne l’est pas).
- Le dual `sqrt(−1)` a `norm = 1`, ce **n’est donc pas** une erreur de
  domaine. L’itération peut produire une partie imaginaire énorme. C’est le
  comportement d’origine.

Guesses initiales : `1` dans le demi-plan droit, `exp(i π / 6)` (ou le
conjugué) sinon pour la racine carrée ; `exp(i π / 9)` pour la cubique.

### Exponentielle

Si l’argument est réel, `exp` de l’hôte. Si `norm_norm(z) > 1`,
`exp(z) = exp(z/2)²` (multiplication de la bibliothèque). Sinon une série
de Horner sur la partie imaginaire, puis multiplication par `exp(real)`.
20 termes par défaut.

Conséquence : `exp(i π) ≈ −1` dans le système ordinaire,
`exp(j) = cosh 1 + j sinh 1` dans le système split, `exp(ε) = 1 + ε` dans
le système dual.

`polar(ρ, θ)` vaut `exp(log ρ + θ · unité)`, avec **cet** `exp`, pas
`std::polar`.

### Logarithme et argument

Si `z = 0` ou `norm(z) < 0`, erreur de domaine. Les réels positifs
utilisent `log` de l’hôte. Sinon on ramène dans une bande par décalages
d’exposant (en ajoutant `ln 2`), puis une identité AGM :

```
log(z) ≈ (π/2) * (1/AGM(1, q) − 1/AGM(1, z q))
```

avec un `q` minuscule. 40 pas externes et 8 pas AGM par défaut. `arg(z)`
est `imag(log(z))`.

### Cosinus et sinus

Cosinus : réduction de la partie réelle modulo `2π`, diminutions de moitié
jusqu’à `norm_norm ≤ 1` (`cos 2α = 2 cos²α − 1`), puis série. 20 termes par
défaut.

**Le sinus est `cos(π/2 − z)`**, comme dans l’original, pas une série
indépendante.

`tan`, `cot`, `sec`, `csc` sont des rapports de ceux-là.

### Fonctions hyperboliques

```
sinh z = (exp(z) − exp(−z)) / 2
cosh z = (exp(z) + exp(−z)) / 2
```

et les rapports habituels. Les inverses hyperboliques sont les identités
log/sqrt usuelles, évaluées avec le `log` et le `sqrt` de cette
bibliothèque.

### Puissance

Les exposants entiers réels passent par l’exponentiation binaire
(`sc_pow_int`). Sinon `exp(log(base) * exposant)`. Les racines complexes
sont `exp(log(base) / exposant)`, ou Newton `root` lorsque l’exposant est
entier.

---

## Fonctions trigonométriques inverses

Il n’y a pas de série `asin` sous forme close. `acos` est l’inversion de
Newton de `cos` (`invert` / `sc_invert`). Les autres sont des identités
bâties dessus :

```
asin(z) = π/2 − acos(z)
atan(z) = acos(1 / sqrt(z² + 1))
acot(z) = π/2 − atan(z)
asec(z) = acos(1/z)
acsc(z) = π/2 − asec(z)
```

**`atan(−1) = +π/4`**, pas `−π/4`. L’identité perd le signe sur l’axe réel
négatif. C’est l’original, pas un bogue de portage.

`invert` différencie la fonction cible deux fois par pas (deux différences
finies) et s’arrête lorsque le `abs` composante par composante de la mise
à jour passe sous `10⁻³`, ou après `max_iterations` (100 par défaut). Si
une différentiation lance, le C++ relance `complex_convergence_error`.

---

## Gamma et bêta

Produit de Weierstrass, **100 000** itérations par défaut (original).
Réflexion pour `Re(z) < 1/2`, récurrence pour `Re(z) > 1`. Les entiers non
positifs sur l’axe réel sont une erreur de domaine.
`beta(a, b) = Γ(a) Γ(b) / Γ(a+b)`.

---

## Calcul différentiel et intégral

Dérivée par différences finies en `a` :

```
d_real = (f(a + δ) − f(a)) / δ
d_imag = (f(a + i δ) − f(a)) / (i δ)
```

Si les deux ne s’accordent pas (`abs` composante par composante de l’écart
relatif `> 10⁻³`), la fonction n’est pas traitée comme
complexe-différentiable : le C écrit `SC_ERR_DOMAIN`, C++/Haskell lancent.
Sinon le résultat est la moyenne. `δ = 10⁻⁶` par défaut.

L’intégrale de `a` à `b` est une somme de Riemann à gauche de pas
`δ (b − a)`, tant que `δ · k ≤ 1`. `δ = 10⁻⁶` par défaut.

Les callbacks C ont le type `sc_func_s` :
`sc_complex (*)(sc_complex, double system)`. C++ et Haskell acceptent tout
foncteur / fonction compatible.

---

## Hyperopérations et M / T / B

`hyper_xexp(z, n)` applique `exp` `n` fois (ou `log` `|n|` fois si
`n < 0`) puis multiplie par `z`.

`hyper_omega(z, n)` est l’inverse en `z` (Newton `invert`) pour `n ≥ 0` ;
un `n` négatif applique des `exp` supplémentaires à `hyper_omega(z, −n)`.

Le C range `n` et `system` courants dans des variables static de fichier
pendant que `sc_hyper_omega` appelle `sc_invert`. C’est une autre raison
pour laquelle l’API C n’est pas réentrante.

`M`, `T`, `B` (et les surcharges C++ par lots sur tableaux) sont les
applications de style Mandelbrot d’origine : `M(a, c, f, L)` vaut
`a − f(a) + c` ou `f(a) + c` selon `L` ; `T` conjugue d’abord ; `B` prend
d’abord le `abs` composante par composante. `sc_square` / le `pow(x, 2)`
par défaut en C++ est le quadratique habituel.

---

## Détails de la surcouche C++

```cpp
#include "scomplex.hpp"
using namespace std::complex::experimental::zhangzl;

complex<double>::set_system(-1.0);
complex<double> z(2.0, 1.0);
auto w = exp(z);
complex<double> i("i");
```

- Espace de noms : `std::complex::experimental::zhangzl`.
- **Ne pas** `#include <complex>` comme on le ferait pour `std::complex` ;
  cet en-tête définit déjà `complex` dans un `std::complex` imbriqué.
- Compiler avec `/utf-8` sous MSVC. Les sources sont en UTF-8 ; sans ce
  drapeau, MSVC peut émettre C4819 et avaler la définition de `_system`
  (LNK2019).
- Lier `scomplex.cpp`. Sans cela, `shift_left` / `shift_right` /
  `round_left` / `medium_mod` restent non résolus.
- `COMPLEX_ALWAYS_INLINE` vaut `__forceinline` sous MSVC. Combiné à
  `inline` sur certains opérateurs, cela donne C4141 ; cela n’affecte pas
  l’édition de liens.
- `sqrt` d’une valeur de norme négative lance `std::domain_error` avec
  `"Square root of non-positive number"`, comme l’original, pas
  `complex_domain_error`.

---

## Détails C

```c
#include "scomplex.h"

sc_complex a = sc_make(2.0, 1.0);
sc_complex b = sc_make(1.0, 3.0);
sc_complex p = sc_mul(a, b, SC_SYSTEM_COMPLEX); /* -1 + 7i */
sc_error_t e = sc_get_error();
```

Compiler et lier `scomplex.c`. C11 suffit. `extern "C"` est fourni pour
les appelants C++ de l’API C.

`sc_get_error` efface le drapeau. Pour inspecter sans consommer, utiliser
`sc_peek_error`.

---

## Détails Haskell

```haskell
import qualified SComplex as SC

let z = SC.fromRealImag 2 1
    w = SC.multiply SC.defaultSystem z (SC.fromRealImag 1 3)
```

La plupart des opérations prennent `system` en premier. `add` /
`subtract` / `conjugate` / `fromString` / `abs` non, car elles sont
indépendantes du système.

`sqrt`, `exp`, `log`, `sin`, … du Prelude sont masqués ; utiliser
`SComplex` ou `qualified Prelude as P`.

---

## Constantes

Fonctions C `sc_pi`, `sc_2pi`, `sc_pi_2`, `sc_pi_3`, `sc_e`, `sc_sqrt2`,
`sc_sqrt1_2`, `sc_cbrt2`, `sc_cbrt1_2`. Statiques C++ sur `complex<_Tp>`
(`pi()`, `e()`, …) plus `i()`. Locaux Haskell `piVal`, `sqrt2`, …

`sc_exp_i_pi_1_6` / `sc_exp_i_pi_1_9` sont les germes de Newton
`exp(i π / 6)` et `exp(i π / 9)`.

---

## Sémantique d’origine à conserver

Si l’on note par rapport aux valeurs principales des manuels, ceci
ressemble à des bogues. C’est la bibliothèque d’origine :

1. `abs` est par composante, pas le module.
2. `sin(z) = cos(π/2 − z)`.
3. `atan(−1) = +π/4`.
4. Le dual `sqrt(−1)` n’est pas une erreur de domaine (`norm = 1`) ; la
   valeur peut être énorme.
5. Le split `sqrt(j)` **est** une erreur de domaine (`norm = −1`).
6. Gamma utilise 100 000 termes de Weierstrass.
7. L’inversion par différences finies s’arrête à `10⁻³`, donc la
   trigonométrie inverse est grossière.

---

## Compilation (MSVC / GHC)

Les sources sont en UTF-8.

```text
cl /utf-8 /std:c11 scomplex.c your.c
cl /utf-8 /std:c++17 /EHsc scomplex.cpp your.cpp
ghc -i. YourMain.hs
```

Équivalents GCC/Clang : compiler `scomplex.c` ou `scomplex.cpp` avec
l’appelant ; le C++ a besoin des exceptions (`-fexceptions`, en général
activé par défaut). L’en-tête C++ utilise des macros d’alignement et de
prefetch MSVC/GCC ; le SIMD n’est pas requis à l’exécution.

---

## Ordre de lecture suggéré dans les sources

1. Multiplication et `norm` — tout le dessin tient dans ce produit.
2. `exp` / `log` / `sqrt` — comment `system` change vraiment l’analyse.
3. `cos` puis `sin` — l’identité de réduction.
4. `invert` puis `acos` / `atan` — pourquoi la trigonométrie inverse a
   cette allure.
5. Chemins d’erreur dans `scomplex.c`, `throw` dans `scomplex.hpp`,
   `throwError` dans `Scomplex.hs`.
