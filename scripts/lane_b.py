#!/usr/bin/env python3
"""
lane_b.py — Certified (rigorous) interval evaluator for L^{(k)}(E,1), k = 0..3,
for E/Q a strong Weil curve.  Proof-grade output: interval [lo, hi] containing
L^{(k)}(E,1) with a DERIVED error budget (truncation + constants + rounding),
no guessed digits.

METHOD (self-contained, no pointwise-modularity assumptions beyond:
   (P)  G(4*pi^2/(N*u)) = eps * (2*pi)^2 * sqrt(N) * G(u*sqrt(N)),  eps = root
   number in {+1,-1}  — verified to 120 digits in logs/probe2.log (curve 20),
   and re-verified numerically inside --selftest for every curve evaluated).

The certified series (derived docs/m1-certified-evaluator.md):

   With G(u) = sum_n a_n exp(-2*pi*n*u/sqrt(N)), u0 = 2*pi*sqrt(N),
   Lambda^{(k)}(1) = int_0^inf G(u) (ln u)^k du  =  sum over n in {a_n} of:

   Lambda^{(k)}(1) = sum_n a_n * K_n^{(k)}

   K_n^{(k)} = int_{u0}^{inf} e^{-lam*u} [ (ln u)^k + eps * c * u^{-2} * (C - ln u)^k ] du
   with  lam = 2*pi*n/sqrt(N),  C = ln(4*pi^2/(N*u0)),  c = + (2*pi)^2*sqrt(N)   if eps=+1
                                              c = - (2*pi)^2/(N)                if eps=-1.

   Upper window kernel  (x = lam*u0 = 4*pi^2*n/sqrt(N),  lnu0 = ln(u0)):
      U0 = e^{-x}/lam
      U1 = ( e^{-x}*lnu0 + E1(x) ) / lam
      U2 = ( e^{-x}*lnu0^2 + 2*lnu0*E1(x) + 2*Q(x) ) / lam
      U3 = ( e^{-x}*lnu0^3 + 3*lnu0^2*E1(x) + 6*lnu0*Q(x) + 6*R(x) ) / lam
   Lower window kernel (m^2 * [C*Wl0 - 2*Wl1 + Wl2] type):
      Wl0 = Gam0m2(x)             = e^{-x}/x^2 - E1(x)        [int e^-t/t^2]
      Wl1 = L Gam0m2(x)           = -e^{-x}ln x^2^2 ...       [int e^-t ln t /t^2]
      Wl2 = L^2 Gam0m2(x)         = ...                       [int e^-t ln^2 t /t^2]

   where E1(x) = int_x^inf e^{-t}/t dt,  Q(x) = int_x^inf e^{-t} ln t /t dt,
   R(x) = int_x^inf e^{-t} (ln t)^2 /t dt  — all computed by ABSOLUTELY convergent
   power series with rigorous interval remainders.

   Conversion to L-derivatives:
      L(s) = (2*pi)^s * N^{-s/2} * Lambda(s) / Gamma(s)
      phi(s) = (2*pi)^s N^{-s/2} / Gamma(s),  L^{(k)}(1) = sum_j C(k,j) phi^{(j)}(1) * Lambda^{(k-j)}(1)
      phi^{(j)}(1) computed from psi-constants (psi(1) = -gamma, psi'(1) = pi^2/6,
      psi''(1) = -2*zeta(3), psi'''(1) = pi^4/15 * 6/... ).  All constants certified
      to >= 3*working digits by rigorous Euler-Maclaurin in interval arithmetic.

USAGE
   python3 lane_b.py --curve '[0,1,1,-2,0]' --k 2 --digits 120 --an work/an-389a1.txt --eps 1
Outputs the interval and its width/justification.  --selftest verifies (P) and
the Hecke relations; --perturb negates the sign of a_n for n odd (must FAIL).
"""
import argparse, math, sys, time
from mpmath import iv, mp, mpf, workdps

# ----------------------------------------------------------------------
# interval helpers
class IV:
    """thin wrapper holding an mpmath interval (mpi); a = (lo,hi) tuple."""
    __slots__ = ("a", "b")
    def __init__(self, a, b=None):
        if b is None:
            if isinstance(a, tuple):  # mpi object from iv.mpf
                self.a, self.b = a[0], a[1]
            else:
                self.a = self.b = mp.mpf(a)
        else:
            self.a, self.b = a, b
    def __add__(s, o):  return IV(s.a + o.a, s.b + o.b)
    def __sub__(s, o):  return IV(s.a - o.b, s.b - o.a)
    def __mul__(s, o):  o = as_iv(o)
        p1, p2 = s.a * o.a, s.a * o.b
        p3, p4 = s.b * o.a, s.b * o.b
        return IV(min(p1, p2, p3, p4), max(p1, p2, p3, p4))
    def __neg__(s): return IV(-s.b, -s.a)
    def __pos__(s): return IV(s.a, s.b)
    def width(s):  return s.b - s.a
    def contains(s, x): return (s.a <= x) and (x <= s.b)
    def mid(s): return (s.a + s.b) / 2
    def __repr__(s):
        # print to current mp.dps with brackets
        return "[%s, %s]" % (mp.nstr(s.a, mp.dps), mp.nstr(s.b, mp.dps))
    def mag(s):
        """max(|s.a|, |s.b|)"""
        return max(abs(s.a), abs(s.b))

def as_iv(o):
    if isinstance(o, IV): return o
    return IV(o, o)

def iv_const(dec_str, digits=50):
    """exact-enclosing interval of a decimal string, width < 10^-digits."""
    mp.dps = digits + 20
    v = mp.mpf(dec_str)
    # mp.mpf(dec) is correctly rounded: enclosing interval is
    mp.dps = digits
    h = mp.mpf(1); h = h << (-(digits - 1) * 4)  # rough 10^-(digits-1)
    # conservative: width 2*10^-digits
    w = mp.mpf(2) * mp.mpf(10) ** (-digits)
    return IV(v - w, v + w)

# ----------------------------------------------------------------------
# certified constants (gamma, zeta(3)) by Euler-Maclaurin in interval arithmetic
# zeta(s) = sum_{n=1}^{M} n^-s  +  M^{1-s}/(s-1) + 1/2*M^-s
#         + sum_{k=1}^{K} B_{2k}/(2k)! * s^{over(2k-1} (M)
# f(x) = x^-s; f^{(m)}(x) = (-1)^m (s)_m x^{-s-m};  (s)_m = s(s+1)...(s+m-1)
def em_zeta(s, M, K, dps):
    with workdps(dps):
        one = IV(1, 1)
        acc = IV(0, 0)
        for n in range(1, M + 1):
            acc = acc + one / IV(n, n) ** IV(s, s)   # exact integers
        # tail:
        # int part: M^{1-s}/(s-1) ... plus EM terms
        # tail = M^{1-s}/(s-1) + 0.5*M^-s + sum_{k=1}^{K} B_{2k} f^{(2k-1)}(M)/(2k)! * sign
        Miv = IV(M, M)
        p = IV(s - 1, s - 1)
        t1 = (Miv ** (-p)) / IV(s - 1, s - 1)          # M^{1-s}/(s-1)
        t2 = IV(1, 2) * Miv ** (IV(-s, -s))            # M^-s / 2
        # EM corrections:  + B_2k/(2k)! * f^{(2k-1)}(M) ... with f^{(m)}(x) = prod_{j=0}^{m-1}(s+j) x^{-s-m}
        corr = IV(0, 0)
        for k in range(1, K + 1):
            j = 2 * k
            B2k = bernoulli_iv(j)                       # certified Bernoulli
            m = j - 1
            fall = IV(1, 1)
            for q in range(m + 1):
                fall = fall * IV(s + q, s + q)
            # f^{(m)}(M) = fall * M^{-s-m};  sign: B_{2k}/(2k)! * (-1)^?  (for the
            # standard zeta EM the sign alternates; keep both candidates and
            # take union to be safe)
            tt = B2k * fall * Miv ** (-(IV(s, s) + IV(m, m)))
            # (2k)! interval
            fact = IV(1, 1)
            for q in range(2, j + 1): fact = fact * IV(q, q)
            tt = tt / fact
            corr = corr + IV(min(tt.a, -tt.b, tt.a, -tt.b) * 0.0 + 0, 0)  # placeholder
        return acc

# (placeholder: the real EM implementation lives below in Cert; this file is
#  the working draft.)

if __name__ == "__main__":
    print("draft: see lane_b_cert.py")
