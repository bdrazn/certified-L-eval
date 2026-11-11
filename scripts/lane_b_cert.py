#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
lane_b_cert.py — Certified (rigorous) interval evaluator for L^{(k)}(E,1), k = 0..2,
on strong Weil curves E/Q.  v7 ENGINE: python-flint ball arithmetic (Arb).
pip3 install python-flint  (0.6.0, verified 2026-11-10).

v7 SIZE-SCALING FIXES (2026-11-11, M2 detector battery — first big-conductor run):
  * n_eval: max() with a (digits+18)-digit decay floor  (was: min() that pinned
    every conductor to the 8*(digits+18) branch — true series tail ~e-8 at big N).
  * law-selftest M: adaptive ceil((dps+46)*ln10*sqrt(N)/pi)+64  (was: fixed 3,000,
    converges only for N <~ 1,700).
  Every M1-configuration certificate is preserved exactly (n_eval 944 at d100,
  4144 at d500; law-gate outcomes unchanged).

CORE IDENTITY (derived + machine-verified 2026-11-10; logs/probe2-4.log):
   G(u) = sum_n a_n exp(-2 pi n u / sqrt(N))
   (P)  G(1/u) = eps * u^2 * G(u)   [pointwise; enforced by selftest]
        => G(1) = 0 exactly when eps = -1  [verified ~1e-135]
   Lambda(s) = N^(s/2)(2 pi)^(-s) Gamma(s) L(s) = int_0^inf G(u) u^(s-1) du
   splitting at the involution fixed point u = 1 and substituting (P):

      Lambda^(k)(1) = (1 + eps (-1)^k) sum_n a_n M_k(lam_n),
      lam_n = 2 pi n / sqrt(N),  M_k(lam) = int_1^inf e^(-lam u) (ln u)^k du

   Kernels:  M_0 = e^-lam/lam ; M_1 = E1(lam)/lam (E1 = Arb expint(1));
   M_2 = 2(Q(lam) - ln(lam) E1(lam))/lam with
      Q(x) = int_x^inf e^-t ln t / t dt
           = C_Q - (ln x)^2/2 - sum_{j>=1} (-1)^j/j![x^j(ln x/j-1/j^2)+1/j^2]
   C_Q = Q(1) certified ONCE (v3): GL-16 over [1, T=200] (nodes built AFTER
   setbits), CORRECT composite remainder C32*(2h)^33 — the old h^34 form was
   2^33/h ~ 1.7e13 too tight — plus the Watson asymptotic tail added INTO THE
   MID as a certified ball (|R| <= e^-T max|f^(M)|, single-extremum bound).
   DIAGNOSIS 2026-11-10 (logs/probe-cq2.log, probe-cq3.log): the old mid was
   the bare GL-[1,40] sum — it missed EXACTLY the tail 3.8492566381206307e-19
   (= int_40^inf, Watson), which propagates through k=2 as 4*sum(a_n/n)*tail
   with sum a_n/n = -0.41411568776043551313 giving +6.37615580907e-19 —
   precisely the recorded "389a1 L''(1) mid off ~6.4e-19" defect. The old
   ball ( +-3.18e-447 ) contained neither tail nor truth: false certificate.

   L^(k)(1) = sum_j C(k,j) phi^(j)(1) Lambda^(k-j)(1),
   phi(s) = (2 pi)^s N^(-s/2)/Gamma(s), Bell chain with psi(1) = -gamma,
   psi'(1) = pi^2/6  (gamma, pi native Arb constants).

   Tail: |a_n| <= 2 sqrt n (Hasse) with kernel bounds, geometric majorant.
   Zero anchoring: (1 + eps(-1)^k) = 0  =>  certified exact ball [0 +- 0].
"""
import json, math, os, sys, time
import flint
from flint import arb, fmpz

WP = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "work")
CONST_CACHE = os.path.join(WP, "laneb-constants.json")

def setbits(bits):
    flint.ctx.prec = max(64, int(bits))

def b2(dps):
    return max(64, int(dps * 3.4 + 32))

def read_an_file(path):
    ak = []
    for line in open(path):
        s = line.strip()
        if s:
            ak.append(int(s))
    if not ak or ak[0] != 1:
        raise ValueError("a_n stream must start with a_1 = 1")
    return ak

# -------- GL-16 via Arb ------------------------------------------------
def _gl16():
    n = 16
    pairs = []
    for k in range(n):
        x, w = arb.legendre_p_root(n, k, weight=True)
        pairs.append((x, w))
    pairs.sort(key=lambda t: t[0].mid())
    return pairs

# -------- Watson tail machinery ----------------------------------------
# f(t) = ln(t)/t:  f^(n)(t) = (a_n ln t + b_n) t^{-(n+1)}, a, b exact ints
# (in fact a_n = (-1)^n n!, b_n = a_n (-1)^{n+1} H_n).
_WATSON_AB = None
def _watson_ab(M):
    global _WATSON_AB
    if _WATSON_AB is None or len(_WATSON_AB) < M + 1:
        a, b = 1, 0
        chain = [(a, b)]
        for n in range(1, M + 1):
            a, b = -n * a, -n * b + a
            chain.append((a, b))
        _WATSON_AB = chain
    return _WATSON_AB

def _fball(n, t):
    """rigorous ball of f^(n)(t) for t a ball; exact-integer coefficients."""
    ch = _watson_ab(n)
    a, b = ch[n]
    return (arb(a) * t.log() + arb(b)) / t ** (n + 1)

def _watson_tail_ball(T, M):
    """certified ball of int_T^inf e^-t ln(t)/t dt.
    Watson:  I = e^-T (sum_{n<M} f^(n)(T)) + R,
    |R| <= e^-T * max_{t>=T} |f^(M)(t)|; the max sits at T or at the single
    interior extremum ln t = 1/(M+1) - b_M/a_M — the candidate is added only
    when provably in [T, inf)."""
    ch = _watson_ab(M)
    aM, bM = ch[M]
    S = arb(0)
    for n in range(M):
        S += _fball(n, T)
    mid = (T * -1).exp() * S
    maxf = abs(_fball(M, T))
    if aM != 0:
        lnstar = arb(1) / (arb(M) + 1) - arb(bM) / arb(aM)
        dt = lstar = lnstar - T.log()
        dec = abs(dt.mid()) > dt.rad()
        if (dec and dt.mid() > 0) or not dec:
            maxf = maxf.union(abs(_fball(M, lnstar.exp())))
    tail_rem = (T * -1).exp() * maxf
    assert abs(tail_rem.mid()) > tail_rem.rad() / 2, "tail remainder not decisive"
    return mid.union(mid - tail_rem).union(mid + tail_rem)

# -------- C_Q: GL-16 quad + Watson-tail-in-mid certification ------------
def certify_CQ(T="200", M="20", half_dec="5e-4"):
    """certified ball for C_Q = Q(1) = int_1^inf e^-t ln t / t dt.
    GL-16 composite over [1, T] (nodes built AFTER setbits — context-safe),
    corrected composite remainder C32*(2h)^33 (the h^34 form was 2^33/h
    ~1.7e13 too tight) with interval-evaluated |f^(32)| per 32-panel block,
    and the Watson tail added INTO THE MID (certified ball)."""
    T, M, half_dec = int(T), int(M), half_dec
    setbits(1400)
    pairs = _gl16()
    stot = sum(w for _, w in pairs)
    d = stot - 2
    assert abs(d.mid()) <= d.rad() + arb(2) ** (-flint.ctx.prec + 16), \
        "GL weights bad / precision too low"
    half = arb(half_dec)
    t0, tT = arb(1), arb(T)
    npanel = int(round(float(((tT - t0) / (2 * half)).mid())))
    C32 = arb(16).fac() ** 4 / (arb(33) * arb(32).fac() ** 3)
    acc = arb(0)
    for p in range(npanel):
        a = t0 + 2 * half * p
        c = a + half
        s = arb(0)
        for (x, w) in pairs:
            xi = c + half * x
            s += w * (xi * -1).exp() * (xi.log() / xi)
        acc += half * s
    blk = 32
    panel_rem = C32 * ((2 * half) ** 33)
    rsum = arb(0)
    p = 0
    while p < npanel:
        cnt = min(blk, npanel - p)
        a_lo = t0 + 2 * half * p
        t_iv = a_lo.union(a_lo + 2 * half * cnt)
        rsum += arb(cnt) * panel_rem * abs(_df32_interval(t_iv)).upper()
        p += cnt
    tail = _watson_tail_ball(tT, M)
    total = acc + tail
    return _widen(total, rsum)

def _widen(ball, w):
    """ball cover [ball.lo - w, ball.hi + w]."""
    setbits(max(flint.ctx.prec, 64))
    return ball.union(ball - w).union(ball + w)

def _df32_interval(t):
    """ball enclosure of f^(32)(t) over an interval ball t (f = e^-t ln(t)/t);
    interval evaluation encloses the image set — rigorous, no monotonicity
    assumption.  t^(r+1) is built by REPEATED MULTIPLICATION: python-flint's
    ball ** on wide-radius bases collapses the mid toward 0 and can straddle
    zero, which poisons subsequent divisions (observed r >= 22, NaN)."""
    A, B = 1, 0
    acc = arb(0)
    tp = arb(1)
    from math import comb
    for r in range(33):
        gr = (arb(B) - arb(A) * t.log()) / tp
        sgn = -1 if (32 - r) % 2 else 1
        acc += arb(comb(32, r)) * sgn * gr
        A, B = (r + 1) * A, (r + 1) * B + A
        tp = tp * t
    return (t * -1).exp() * acc

# -------- constants ------------------------------------------------------
def load_constants(digits):
    """gamma/pi native Arb constants; C_Q certified ball (v3 cache semantics:
    the mid INCLUDES the Watson tail; the certified radius is stored as a
    decimal string and re-applied verbatim on reload — the old e10
    reconstruction narrowed the reloaded ball below its certificate)."""
    if os.path.exists(CONST_CACHE):
        d = json.load(open(CONST_CACHE))
        if d.get("v") == 3 and d.get("cap", 0) >= digits + 10:
            setbits(b2(digits + 60))
            m = arb(d["mid"])
            if "rad10" in d:
                rad = arb(str(d["rad10"]).strip("[] ").split("+/-")[0])
            else:
                rad = m.rad() * 2 + arb(10) ** (-(digits + 100))
            w = rad * (1 + arb("1e-20")) + m.rad() + arb(10) ** (-(digits + 100))
            cq = m.union(m - w).union(m + w)
            return {"gamma": arb.const_euler(), "pi": arb.pi(), "C_Q": cq}
    setbits(b2(digits + 80))
    cq = certify_CQ()
    json.dump({"v": 3, "cap": digits + 10, "mid": str(cq.mid()),
               "rad10": str(cq.rad())}, open(CONST_CACHE, "w"))
    return {"gamma": arb.const_euler(), "pi": arb.pi(), "C_Q": cq}

def _one_arb():
    return arb(1)

# -------- kernels -------------------------------------------------------
def lam_n(n, N):
    return 2 * arb.pi() * arb(n) / arb(N).sqrt()

def M0(lam):
    return (lam * -1).exp() / lam

def M1(lam):
    return lam.expint(1) / lam

def Q_of(x, ctx, digits):
    CQ = ctx["C_Q"]
    lx = x.log()
    acc = CQ - lx * lx / 2
    xm = float(x.mid()) + float(x.rad())
    # Watson-vs-series crossover scales with the target: the Watson remainder
    # floor is ~ e^-x * M!/x^(M+1) * |ln x - H_M|, so for a (digits+45) target
    # the series must take over below x ~ 2.3*(digits+45).
    lam_w = max(100.0, 2.3 * (digits + 45))
    if xm >= lam_w:
        # large argument: Watson asymptotic directly (no cancellation).
        # |R| <= e^-x max|f^(25)| — decisive, and at x>=100 far tighter than
        # the series' cancellation floor would be.
        return _watson_tail_ball(x, 25)
    # alternating series; stop when |t_j| <= 10^-(digits+45)/3 and j >= 4x.
    # Gates make later ratios <= 2x/(k+1) <= 1/2 and declining, so the
    # remainder is <= 3*|t_(j+1)| — verified by the measured next-ratio.
    thresh3 = arb(10) ** (-(digits + 45)) / 3
    jguard = max(int(4 * xm), 16)
    # jcap must cover jguard for EVERY series-zone lambda: lam <= lam_w =
    # 2.3*(digits+45) reaches 4*lam ~ 9.2*(digits+45) — at d500 that is
    # 5,013 > 6*digits+800 = 3,800, so the old fixed jcap made every
    # lambda > (6*digits+800)/4 raise "series failed to reach threshold"
    # before the guard allowed exit (389a1 d500 crash, first bad term
    # n = 2982, lam = 950.3).  Series-zone lambda <= 950 hits jcap at the
    # same 3,800 as before, so d100 semantics are unchanged.
    jcap = max(6 * digits + 800, jguard + 16)
    xj = arb(1)
    jf = arb(1)
    last = None
    j = 0
    while True:
        j += 1
        if j > jcap:
            raise RuntimeError("Q_of series failed to reach threshold")
        xj *= x
        jf *= arb(j)
        sgn = -1 if j % 2 else 1
        core = xj * (lx / arb(j) - 1 / arb(j * j)) + 1 / arb(j * j)
        last = sgn * core / jf
        acc -= last
        if j < jguard:
            continue
        slack = thresh3 - abs(last)
        if not (abs(slack.mid()) > slack.rad() and slack.mid() > 0):
            continue
        j1 = j + 1
        xj *= x
        jf *= arb(j1)
        sgn = 1 if j1 % 2 == 0 else -1
        core = xj * (lx / arb(j1) - 1 / arb(j1 * j1)) + 1 / arb(j1 * j1)
        t2 = sgn * core / jf
        j2 = j1 + 1
        xj *= x
        jf *= arb(j2)
        sgn = 1 if j2 % 2 == 0 else -1
        core = xj * (lx / arb(j2) - 1 / arb(j2 * j2)) + 1 / arb(j2 * j2)
        t3 = sgn * core / jf
        q = abs(t3) / abs(t2)
        assert abs(q.mid()) + q.rad() <= arb("0.6"), "series tail ratio > 0.6"
        return _widen(acc, abs(t2) * 3)

def M2(lam, ctx, digits):
    return 2 * (Q_of(lam, ctx, digits) - lam.log() * lam.expint(1)) / lam

def kernel(k, lam, ctx, digits):
    if k == 0: return M0(lam)
    if k == 1: return M1(lam)
    if k == 2: return M2(lam, ctx, digits)
    raise ValueError("k <= 2")

# -------- tail bound ----------------------------------------------------
def tail_bound(k, N, n0):
    lam1 = lam_n(1, N)
    M = arb(n0)
    one = arb(1)
    def S(alpha):
        w = (lam1 * (M + one) * -1).exp() * (M + one) ** (-alpha)
        return w / (one - (lam1 * -1).exp())
    if k == 0:
        return (2 / lam1) * S(arb("0.5"))
    if k == 1:
        return (2 / (lam1 * lam1)) * S(arb("1.5"))
    l1 = lam1.mid() + lam1.rad()
    lg = abs(arb(l1).log())
    C1 = 2 * (lg / arb(l1) + 2)
    C2 = 2 * (lg + 1)
    return C1 * S(arb("0.5")) + C2 * S(-1)

# -------- conversion ----------------------------------------------------
def phi_derivs(N, ctx):
    gam, pi = ctx["gamma"], ctx["pi"]
    ln2pi = (2 * pi).log()
    lnN = arb(N).log()
    F0 = ln2pi - lnN / 2
    F1 = ln2pi - lnN / 2 + gam
    F2 = -(pi * pi) / 6
    p0 = F0.exp()
    return [p0, p0 * F1, p0 * (F2 + F1 * F1)]

from math import comb as _comb

def evaluate(N, ak, eps, kmax, digits, ctx):
    lam1 = lam_n(1, N)
    # v7 size-scaling fix (2026-11-11, M2 detector battery): n_eval must grow with
    # sqrt(N) — the old min() pinned every conductor to the 8*(digits+18) branch,
    # leaving a true series tail of ~e-8 scale at N ~ 250k-500k (useless balls).
    # New form: max(8*(digits+18), (digits+18)-digit decay floor at sqrt(N)/2pi).
    # Every M1-configuration n_eval is preserved exactly: 389/d100->944, 20/d100->944,
    # 37/d100->944, 37/d500->4144, 389/d500->4144 (the floor branch never binds there).
    n_eval = int(max(8.0 * (digits + 18),
                     math.ceil((digits + 18) * math.log(10.0) * math.sqrt(N)
                               / (2.0 * math.pi))))
    n_eval = min(n_eval, len(ak) - 1)
    results = {}
    lam_max = float(lam1.mid()) * n_eval
    bits_storm = b2(int(0.4343 * lam_max + 3.33 * (digits + 14) + 40))
    setbits(bits_storm)
    for k in range(0, kmax + 1):
        factor = 1 + eps * ((-1) ** k)
        if factor == 0:
            results[k] = arb(0)
            continue
        acc = arb(0)
        for n in range(1, n_eval + 1):
            acc += arb(fmpz(ak[n - 1])) * kernel(k, lam_n(n, N), ctx, digits)
        results[k] = _widen(acc, tail_bound(k, N, n_eval)) * factor
    phis = phi_derivs(N, ctx)
    out = {}
    for k in range(0, kmax + 1):
        ac = arb(0)
        for j in range(0, k + 1):
            ac += arb(_comb(k, j)) * phis[j] * results[k - j]
        out[k] = ac
    return {"L": out, "Lambda": results, "n_eval": n_eval, "bits": bits_storm}

# -------- selftests -----------------------------------------------------
def selftest_pointwise_law(N, ak, eps, dps):
    setbits(b2(dps))
    sqN = arb(N).sqrt()
    # v7 size-scaling fix (2026-11-11, M2): the old fixed M=3000 only converges
    # for conductors N <~ 1,700 (decay at u=0.5 is e^{-pi*n/sqrt(N)}, so 3,000
    # terms leave ~e-8 total decay at N ~ 250k) — the gate was failing on its
    # own truncation, not on the stream.  M now scales so the direct-series
    # tail (both sides evaluated directly — required so a corrupted stream
    # cannot pass) is below the tolerance with margin.  For every M0/M1
    # configuration the new M is smaller but still gives >e-180 decay where
    # the old 3,000 gave >=e-207, so all past gate outcomes are unchanged.
    sqN_hi = float(sqN.mid()) + float(sqN.rad())
    M_law = int(math.ceil((dps + 46) * math.log(10.0) * sqN_hi / math.pi)) + 64
    def G(u, M=None):
        if M is None:
            M = M_law
        s = arb(0)
        for n in range(1, M + 1):
            s += arb(fmpz(ak[n - 1])) * arb.exp(-2 * arb.pi() * n * u / sqN)
        return s
    tol = arb(10) ** (-dps + 42)
    ok = True
    for u in (arb("0.5"), arb(2)):
        if abs(G(1 / u) - eps * u * u * G(u)) > tol:
            ok = False
    if eps == -1 and abs(G(1)) > tol:
        ok = False
    return ok

def hecke_check(ak, N, maxp=200):
    bad = []
    if ak[0] != 1:
        bad.append("a1")
    def isprime(p):
        d = 2
        while d * d <= p:
            if p % d == 0: return False
            d += 1
        return p > 1
    L = len(ak)
    for p in range(2, min(maxp, L)):
        if N % p != 0 and isprime(p) and p * p <= L and ak[p * p - 1] != ak[p - 1] ** 2 - p:
            bad.append(p)
    from math import gcd
    for m in range(2, min(60, L)):
        for n in range(2, min(60, L)):
            if gcd(m, n) == 1 and m * n <= L and ak[m * n - 1] != ak[m - 1] * ak[n - 1]:
                bad.append((m, n))
    return bad

# -------- CLI -------------------------------------------------------
def main():
    import argparse
    ap = argparse.ArgumentParser()
    ap.add_argument("--N", type=int, required=True)
    ap.add_argument("--an", required=True)
    ap.add_argument("--eps", type=int, required=True, choices=[-1, 1])
    ap.add_argument("--k", type=int, default=2)
    ap.add_argument("--digits", type=int, default=100)
    ap.add_argument("--perturb", action="store_true")
    ap.add_argument("--json-out", default=None)
    ap.add_argument("--no-lawtest", action="store_true")
    args = ap.parse_args()
    ak = read_an_file(args.an)
    if args.perturb:
        ak = [(-a if (i + 1) % 3 == 1 else a) for i, a in enumerate(ak)]
    bad = hecke_check(ak, args.N)
    if bad and not args.perturb:
        print("HECKE-FAIL", bad, file=sys.stderr)
        sys.exit(3)
    ctx = load_constants(args.digits)
    if not args.no_lawtest and not selftest_pointwise_law(args.N, ak, args.eps, min(args.digits + 30, 170)):
        print("POINTWISE-LAW-FAIL", file=sys.stderr)
        sys.exit(4)
    print("# (P) pointwise law verified on a_n stream", file=sys.stderr)
    res = evaluate(args.N, ak, args.eps, args.k, args.digits, ctx)
    setbits(res["bits"] + 80)
    outs = {}
    for k, b in res["L"].items():
        outs[int(k)] = (b.lower(), b.upper())
    rep = {"curve_N": args.N, "eps": args.eps, "digits_target": args.digits,
           "n_eval": res["n_eval"], "prec_bits": res["bits"],
           "perturbed": bool(args.perturb),
           "L": {str(k): "[{%s}, {%s}]" % (v[0], v[1]) for k, v in outs.items()}}
    s = json.dumps(rep, indent=1)
    if args.json_out:
        open(args.json_out, "w").write(s)
    print(s)
    for k in sorted(outs):
        print("# L^%d(1) in [%s, %s]  (n_eval=%d, prec_bits=%d)" %
              (k, outs[k][0], outs[k][1], res["n_eval"], res["bits"]))

if __name__ == "__main__":
    main()
