{-# LANGUAGE DeriveAnyClass #-}
-- | Complex / split-complex / dual-number arithmetic (Haskell).
--
-- 'system' is an explicit argument: -1 ordinary complex (i^2 = -1),
-- +1 split / hyperbolic (j^2 = +1), 0 dual (epsilon^2 = 0).
-- Algebra matches scomplex.hpp / scomplex.h. Errors are thrown as
-- 'ComplexError', not swallowed. Transcendentals use system-aware
-- multiplication (series, Newton, AGM), not Prelude polar formulas.
module SComplex (
    Complex(..),
    ComplexError(..),
    fromRealImag,
    fromString,
    add, subtract, multiply, divide, modulo,
    conjugate, norm, normNorm, abs,
    sqrt, cbrt, root, exp, log, pow,
    sin, cos, tan, cot, sec, csc,
    sinh, cosh, tanh, coth, sech, csch,
    asin, acos, atan, acot, asec, acsc,
    asinh, acosh, atanh, acoth, asech, acsch,
    gamma, beta,
    polar, arg, phase,
    logBaseC,
    differential, invert, integral,
    hyperXexp, hyperOmega,
    mFunc, tFunc, bFunc,
    withSystem,
    defaultSystem
) where

import Control.Exception (Exception, throw)
import Data.Char (isSpace)
import qualified Prelude as P
import Prelude hiding (
    sqrt, exp, log, sin, cos, tan, abs, subtract,
    sinh, cosh, tanh, asin, acos, atan, asinh, acosh, atanh)

data Complex = Complex { real :: !Double, imag :: !Double }
    deriving (Eq, Show)

data ComplexError
    = DivisionByZero
    | Overflow
    | DomainError String
    | ConvergenceError
    | InvalidFormat String
    deriving (Show, Eq, Exception)

defaultSystem :: Double
defaultSystem = -1

withSystem :: Double -> (Double -> a) -> a
withSystem s f = f s

throwError :: ComplexError -> a
throwError = throw

isFinite :: Double -> Bool
isFinite x = not (isNaN x || isInfinite x)

checkOverflow :: Double -> Double
checkOverflow x
    | not (isFinite x) = throwError Overflow
    | P.abs x > 8.988465674311579e307 = throwError Overflow
    | otherwise = x

checkDivisionByZero :: Double -> Double
checkDivisionByZero x
    | P.abs x < 2.220446049250313e-16 = throwError DivisionByZero
    | otherwise = x

fromRealImag :: Double -> Double -> Complex
fromRealImag r i = Complex (checkOverflow r) (checkOverflow i)

zero, one, iUnit :: Complex
zero = fromRealImag 0 0
one  = fromRealImag 1 0
iUnit = fromRealImag 0 1

piVal, pi2, pi1_2, pi1_3, sqrt2, sqrt1_2, cbrt2, cbrt1_2, ln2 :: Double
piVal   = 3.14159265358979323846
pi2     = 6.28318530717958647692
pi1_2   = 1.57079632679489661923
pi1_3   = 1.04719755119659774615
sqrt2   = 1.41421356237309504880
sqrt1_2 = 0.70710678118654752440
cbrt2   = 1.25992104989487316476
cbrt1_2 = 0.79370052598409973737
ln2     = 0.69314718055994530942

expIPi16, expIPi19 :: Complex
expIPi16 = fromRealImag 0.5 0.86602540378443864676
expIPi19 = fromRealImag 0.93969262078590838405 0.34202014332566873304

stripSpace :: String -> String
stripSpace = filter (not . isSpace)

fromString :: String -> Complex
fromString str = parse (stripSpace str)
  where
    parse "i"  = iUnit
    parse "+i" = iUnit
    parse "-i" = fromRealImag 0 (-1)
    parse s =
        case spanNum s of
            (Just r, "") -> fromRealImag r 0
            (Just r, rest@(c:_)) | c == '+' || c == '-' ->
                let sign = if c == '-' then -1 else 1
                    rest' = drop 1 rest
                in case rest' of
                    "i" -> fromRealImag r (sign * 1)
                    _ -> case spanNum rest' of
                        (Just im, "i") -> fromRealImag r (sign * im)
                        _ -> throwError (InvalidFormat str)
            _ -> throwError (InvalidFormat str)

    spanNum :: String -> (Maybe Double, String)
    spanNum s = case reads s of
        [(x, more)] -> (Just x, more)
        _ -> (Nothing, s)

add :: Complex -> Complex -> Complex
add (Complex r1 i1) (Complex r2 i2) = fromRealImag (r1 + r2) (i1 + i2)

subtract :: Complex -> Complex -> Complex
subtract (Complex r1 i1) (Complex r2 i2) = fromRealImag (r1 - r2) (i1 - i2)

multiply :: Double -> Complex -> Complex -> Complex
multiply system (Complex r1 i1) (Complex r2 i2) =
    fromRealImag (r1 * r2 + system * i1 * i2) (r1 * i2 + i1 * r2)

divide :: Double -> Complex -> Complex -> Complex
divide system (Complex r1 i1) (Complex r2 i2) =
    let denom = checkDivisionByZero (r2 * r2 - system * i2 * i2)
    in fromRealImag ((r1 * r2 - system * i1 * i2) / denom)
                    ((i1 * r2 - r1 * i2) / denom)

roundLeft :: Double -> Double
roundLeft a = a - fromInteger (round a)

modulo :: Double -> Complex -> Complex -> Complex
modulo system a b =
    let Complex qr qi = divide system a b
        q = fromRealImag (roundLeft qr) (roundLeft qi)
    in multiply system q b

conjugate :: Complex -> Complex
conjugate (Complex r i) = Complex r (-i)

norm :: Double -> Complex -> Double
norm system (Complex r i) = r * r - system * i * i

normNorm :: Complex -> Double
normNorm (Complex r i) = r * r + i * i

abs :: Complex -> Complex
abs (Complex r i) = Complex (P.abs r) (P.abs i)

shiftLeft, shiftRight :: Complex -> Integer -> Complex
shiftLeft  (Complex r i) n = fromRealImag (scaleFloat (fromIntegral n) r) (scaleFloat (fromIntegral n) i)
shiftRight z n = shiftLeft z (-n)

powInt :: Double -> Complex -> Int -> Complex
powInt _ _ 0 = one
powInt _ z _ | real z == 0 && imag z == 0 = zero
powInt system z n
    | n < 0 = divide system one (powInt system z (-n))
    | odd n = multiply system z (powInt system z (n - 1))
    | otherwise = let h = powInt system z (n `div` 2) in multiply system h h

sqrtInit :: Double -> Complex -> Complex -> Int -> Complex
sqrtInit system a b0 maxIt = go b0 b0 0
  where
    go b bp k
        | k >= maxIt = b
        | otherwise =
            let nx = divide system (add bp (divide system a bp)) (fromRealImag 2 0)
                Complex dr di = abs (subtract nx bp)
                b' = nx
            in if dr < 2.220446049250313e-16 && di < 2.220446049250313e-16
               then b'
               else go b' b' (k + 1)

sqrt :: Double -> Complex -> Complex
sqrt system a
    | imag a == 0 && real a >= 0 = fromRealImag (P.sqrt (real a)) 0
    | norm system a < 0 = throwError (DomainError "Square root of non-positive number")
    | normNorm a > 2 = multiply system (sqrt system (shiftRight a 1)) (fromRealImag sqrt2 0)
    | normNorm a < 0.5 = multiply system (sqrt system (shiftLeft a 1)) (fromRealImag sqrt1_2 0)
    | real a > 0 = sqrtInit system a one 100
    | imag a > 0 = sqrtInit system a expIPi16 100
    | otherwise = sqrtInit system a (conjugate expIPi16) 100

cbrtInit :: Double -> Complex -> Complex -> Int -> Complex
cbrtInit system a b0 maxIt = go b0 b0 0
  where
    go b bp k
        | k >= maxIt = b
        | otherwise =
            let nx = divide system (add (multiply system bp (fromRealImag 2 0))
                                        (divide system a (multiply system bp bp)))
                                   (fromRealImag 3 0)
                Complex dr di = abs (subtract nx bp)
            in if dr < 2.220446049250313e-16 && di < 2.220446049250313e-16
               then nx
               else go nx nx (k + 1)

cbrtReal :: Double -> Double
cbrtReal x = signum x * (P.abs x ** (1/3))

cbrt :: Double -> Complex -> Complex
cbrt system a
    | imag a == 0 = fromRealImag (cbrtReal (real a)) 0
    | normNorm a > 2 = multiply system (cbrt system (shiftRight a 1)) (fromRealImag cbrt2 0)
    | normNorm a < 0.5 = multiply system (cbrt system (shiftLeft a 1)) (fromRealImag cbrt1_2 0)
    | real a > 0 = cbrtInit system a one 100
    | imag a > 0 = cbrtInit system a expIPi19 100
    | otherwise = cbrtInit system a (conjugate expIPi19) 100

root :: Double -> Complex -> Int -> Complex
root _ a 1 = a
root system a 2 = sqrt system a
root system a 3 = cbrt system a
root system a n
    | imag a == 0 = fromRealImag (real a ** (1 / fromIntegral n)) 0
    | normNorm a > 2 = multiply system (root system (shiftRight a 1) n)
                                       (fromRealImag (2 ** (1 / fromIntegral n)) 0)
    | normNorm a < 0.5 = multiply system (root system (shiftLeft a 1) n)
                                         (fromRealImag (0.5 ** (1 / fromIntegral n)) 0)
    | otherwise = go initG initG (0 :: Int)
  where
    ang = pi1_3 / fromIntegral n
    initG
        | real a > 0 = one
        | imag a > 0 = fromRealImag (P.cos ang) (P.sin ang)
        | otherwise = fromRealImag (P.cos ang) (-P.sin ang)
    go b bp k
        | k >= (100 :: Int) = b
        | otherwise =
            let nx = divide system (add (multiply system bp (fromRealImag (fromIntegral (n - 1)) 0))
                                        (divide system a (powInt system bp (n - 1))))
                                   (fromRealImag (fromIntegral n) 0)
                Complex dr di = abs (subtract nx bp)
            in if dr < 2.220446049250313e-16 && di < 2.220446049250313e-16
               then nx
               else go nx nx (k + 1)

expIter :: Double -> Complex -> Int -> Complex
expIter system a maxIt
    | imag a == 0 = fromRealImag (P.exp (real a)) 0
    | normNorm a > 1 =
        let h = expIter system (shiftRight a 1) maxIt
        in multiply system h h
    | otherwise =
        let temp = fromRealImag 0 (imag a)
            b = foldl (\acc k -> add (divide system (multiply system temp acc)
                                                    (fromRealImag (fromIntegral (k + 1)) 0))
                                     one)
                      one [maxIt, maxIt-1 .. 0]
        in multiply system b (fromRealImag (P.exp (real a)) 0)

exp :: Double -> Complex -> Complex
exp system z = expIter system z 20

agm :: Double -> Complex -> Complex -> Int -> Complex
agm system a0 b0 maxIt = go a0 b0 (add a0 b0) (multiply system a0 b0) 0
  where
    go a _ s p k
        | k >= maxIt = a
        | otherwise =
            let a' = divide system s (fromRealImag 2 0)
                b' = sqrt system p
            in go a' b' (add a' b') (multiply system a' b') (k + 1)

logIter :: Double -> Complex -> Int -> Complex
logIter system a maxIt
    | real a == 0 && imag a == 0 = throwError (DomainError "Logarithm of zero is undefined")
    | imag a == 0 && real a > 0 = fromRealImag (P.log (real a)) 0
    | norm system a < 0 = throwError (DomainError "Logarithm of non-positive number")
    | normNorm a > 1 = add (logIter system (shiftRight a 1) maxIt) (fromRealImag ln2 0)
    | normNorm a < 0.25 = subtract (logIter system (shiftLeft a 1) maxIt) (fromRealImag ln2 0)
    | otherwise =
        let q = iterate (\x -> divide system x (fromRealImag 2 0)) one !! (maxIt - 1)
            ans1 = agm system one q 8
            ans2 = agm system one (multiply system a q) 8
        in if normNorm ans1 < 2.220446049250313e-16 || normNorm ans2 < 2.220446049250313e-16
           then throwError ConvergenceError
           else multiply system (subtract (divide system one ans1) (divide system one ans2))
                                (fromRealImag pi1_2 0)

log :: Double -> Complex -> Complex
log system z = logIter system z 40

logBaseC :: Double -> Complex -> Complex -> Complex
logBaseC system a b = divide system (log system a) (log system b)

arg :: Double -> Complex -> Double
arg system a = imag (log system a)

phase :: Double -> Complex -> Double
phase = arg

polar :: Double -> Double -> Double -> Complex
polar system rho theta = exp system (fromRealImag (P.log rho) theta)

cosIter :: Double -> Complex -> Int -> Complex
cosIter system a0 maxIt
    | imag a0 == 0 = fromRealImag (P.cos (real a0)) 0
    | otherwise =
        let a = reduce a0
        in if normNorm a > 1
           then let h = cosIter system (shiftRight a 1) maxIt
                in subtract (multiply system (multiply system h h) (fromRealImag 2 0)) one
           else foldl (\b i -> step a i b) one [maxIt, maxIt-1 .. 0]
  where
    reduce z
        | real z > pi2 = reduce (fromRealImag (real z - pi2) (imag z))
        | real z < -pi2 = reduce (fromRealImag (real z + pi2) (imag z))
        | otherwise = z
    step a i b =
        let p = divide system (multiply system a b) (fromRealImag (fromIntegral (i + 1)) 0)
        in if odd i then p
           else if i `mod` 4 /= 0 then subtract p one
           else add p one

cos :: Double -> Complex -> Complex
cos system z = cosIter system z 20

sin :: Double -> Complex -> Complex
sin system a = cos system (subtract (fromRealImag pi1_2 0) a)

tan, cot, sec, csc :: Double -> Complex -> Complex
tan system a = divide system (sin system a) (cos system a)
cot system a = divide system (cos system a) (sin system a)
sec system a = divide system one (cos system a)
csc system a = divide system one (sin system a)

sinh, cosh, tanh, coth, sech, csch :: Double -> Complex -> Complex
sinh system a = divide system (subtract (exp system a) (exp system (fromRealImag (-real a) (-imag a)))) (fromRealImag 2 0)
cosh system a = divide system (add (exp system a) (exp system (fromRealImag (-real a) (-imag a)))) (fromRealImag 2 0)
tanh system a = divide system (sinh system a) (cosh system a)
coth system a = divide system (cosh system a) (sinh system a)
sech system a = divide system one (cosh system a)
csch system a = divide system one (sinh system a)

pow :: Double -> Complex -> Complex -> Complex
pow _ _ e | real e == 0 && imag e == 0 = one
pow _ b _ | real b == 0 && imag b == 0 = zero
pow system base expo
    | imag expo == 0 && real expo == fromInteger (round (real expo))
        = powInt system base (round (real expo))
    | otherwise = exp system (multiply system (log system base) expo)

asinh, acosh, atanh, acoth, asech, acsch :: Double -> Complex -> Complex
asinh system a = log system (add a (sqrt system (add (multiply system a a) one)))
acosh system a = log system (add a (sqrt system (subtract (multiply system a a) one)))
atanh system a = divide system (log system (divide system (add one a) (subtract one a))) (fromRealImag 2 0)
acoth system a = divide system (log system (divide system (add a one) (subtract a one))) (fromRealImag 2 0)
asech system a = log system (divide system (add (sqrt system (subtract one (multiply system a a))) one) a)
acsch system a = log system (divide system (add (sqrt system (add one (multiply system a a))) one) a)

differential :: Double -> Complex -> (Complex -> Complex) -> Double -> Complex
differential system a func delta =
    let dReal = divide system (subtract (func (add a (fromRealImag delta 0))) (func a)) (fromRealImag delta 0)
        dImag = divide system (subtract (func (add a (fromRealImag 0 delta))) (func a)) (fromRealImag 0 delta)
        d = abs (divide system (subtract dReal dImag) (func a))
    in if real d > 1e-3 || imag d > 1e-3
       then throwError (DomainError "The function is not differentiable")
       else divide system (add dReal dImag) (fromRealImag 2 0)

invert :: Double -> Complex -> (Complex -> Complex) -> Int -> Complex
invert system a func maxIt = go one one 0
  where
    go b bp k
        | k >= maxIt = b
        | otherwise =
            let d1 = differential system bp func 1e-6
                d2 = differential system bp func 1e-6
                b' = add (subtract bp (divide system (func bp) d1)) (divide system a d2)
                Complex dr di = abs (subtract b' bp)
            in if dr < 1e-3 && di < 1e-3 then b' else go b' b' (k + 1)

acos, asin, atan, acot, asec, acsc :: Double -> Complex -> Complex
acos system a = invert system a (cos system) 100
asin system a = subtract (fromRealImag pi1_2 0) (acos system a)
atan system a = acos system (divide system one (sqrt system (add (powInt system a 2) one)))
acot system a = subtract (fromRealImag pi1_2 0) (atan system a)
asec system a = acos system (divide system one a)
acsc system a = subtract (fromRealImag pi1_2 0) (asec system a)

integral :: Double -> Complex -> Complex -> (Complex -> Complex) -> Double -> Complex
integral system a0 b func delta = go a0 zero (multiply system (func a0) d) (0 :: Int)
  where
    d = multiply system (subtract b a0) (fromRealImag delta 0)
    go a ans ins k
        | delta * fromIntegral k > 1 = ans
        | otherwise =
            let a' = add a d
            in go a' (add ans ins) (multiply system (func a') d) (k + 1)

gammaIter :: Double -> Complex -> Int -> Complex
gammaIter _ a _ | a == one = one
gammaIter system a maxIt
    | imag a == 0 && real a < 0 && fromInteger (round (real a)) == real a
        = throwError (DomainError "Gamma of non-positive integer is undefined")
    | real a < 0.5 =
        divide system (fromRealImag piVal 0)
            (multiply system (sin system (multiply system a (fromRealImag piVal 0)))
                             (gammaIter system (subtract one a) maxIt))
    | real a > 1 =
        multiply system (gammaIter system (subtract a one) maxIt) (subtract a one)
    | otherwise =
        let bStart = let l = log system a in fromRealImag (-real l) (-imag l)
            step i b =
                let b1 = subtract b (log system (add one (divide system a (fromRealImag (fromIntegral i) 0))))
                in add b1 (multiply system a (fromRealImag (P.log (1 + 1 / fromIntegral i)) 0))
        in exp system (foldl (flip step) bStart [1..maxIt])

gamma :: Double -> Complex -> Complex
gamma system a = gammaIter system a 100000

beta :: Double -> Complex -> Complex -> Complex
beta system a b = divide system (multiply system (gamma system a) (gamma system b)) (gamma system (add a b))

hyperXexp :: Double -> Complex -> Int -> Complex
hyperXexp system a n =
    let pos = foldl (\acc _ -> exp system acc) a [1..n]
        neg = foldl (\acc _ -> log system acc) a [1..(-n)]
        ans = if n >= 0 then pos else neg
    in multiply system ans a

hyperOmega :: Double -> Complex -> Int -> Complex
hyperOmega system a n
    | n >= 0 = invert system a (\x -> hyperXexp system x n) 100
    | otherwise =
        let ans = hyperOmega system a (-n)
        in foldl (\acc _ -> exp system acc) ans [1..(-n)]

mFunc :: Double -> Complex -> Complex -> (Complex -> Complex) -> Bool -> Complex
mFunc _ a c func l = add (if l then subtract a (func a) else func a) c

tFunc :: Double -> Complex -> Complex -> (Complex -> Complex) -> Bool -> Complex
tFunc system a c func l = mFunc system (conjugate a) c func l

bFunc :: Double -> Complex -> Complex -> (Complex -> Complex) -> Bool -> Complex
bFunc system a c func l = mFunc system (abs a) c func l
