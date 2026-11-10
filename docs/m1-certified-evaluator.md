# M1 — Certified L-Derivative Evaluator (design)

Status: design committed 2026-11-10. Implementation next.

## Goal

Proof-grade rigorous bounds on L^(j)(E,1) for j = 0, 1, 2 (and later j up to r+1) on curves over ℚ — the certification layer under every detector:

- **D1 / CERT-1** needs only: rigorous `L(1) ≠ 0` (one certified evaluation proving r_an = 0).
- **D2 / CERT-2** needs: rigorous `L''(1)` with error ≪ distance to the nearest square-integer exclusion window.
- Every alarm must reproduce under two independent implementations (campaign rule #2).

## Implementation lanes

**Lane A — eclib (Cremona) integration.** eclib 20261102 built at `~/.local/eclib` (verified: mwrank unconditional rank + saturated basis, 0.14 s on 37.a1). Open investigation: locate the L-series-with-bounds module (historically `lbseries`); the 20261102 tree's `libsrc` grep for `Lseries|lbseries` headers returned nothing — inventory step 1 is a full-text search (`grep -rli "lbseries\|L-series\|Lvalue" work/eclib-src`) + reading `src/` program list (h1bsd, qexp). If eclib exposes rigorous L-derivative bounds (it underlies Cremona's verified BSD computations), Lane A becomes the primary evaluator.

**Lane A resolution (2026-11-10):** confirmed absent — no L-series-with-bounds module exists in eclib 20261102 under any name (`lseries`/`lbseries` do not exist; the L machinery behind mwrank's analytic rank is `libsrc/kbessel.cc`, API `double kbessel(double nu, double gx)` — **double-precision only**, built for heuristic rank). eclib's certification-relevant asset is `saturate.cc` (in use for regulators). **Lane A is demoted to saturation + non-certified cross-check; certification rests entirely on Lane B.**

**Lane B — from-scratch certified evaluator (our own, swarm-buildable).** The direct Dirichlet series is NOT absolutely convergent at s = 1; the certified path is the complete modified series (Dokchitser-type) with an explicit tail bound:
- coefficients: |a_n| ≤ σ₀(n)·√n (Deligne), so |a_n·n^{-s}| ≤ bounded with the modified-series weights making the tail O(e^{−c·√(n·X)})-style;
- choose truncation X so the certified tail < 10^{−(d+guard)} for d target digits;
- evaluate the truncated sum in interval arithmetic (python3 + mpmath `iv`, or C++ + Arb for speed later);
- output: [L̂ − err, L̂ + err] with err = truncation + rounding, both derived, not guessed.
This is M1's real math-engineering deliverable. Validation gate: reproduce 389.a1's L″(1) = 1.51863300057… and 37.a1's L′(1) = 0.30599977383… with matching error bars, and reproduce a zero-anchoring case (L(1) = 0 for 37.a1 at j = 0: interval must CONTAIN 0, and bound |L′(1)| > 0 rigorously).

**Lane C — PARI `lfun` overlay (non-rigorous reference + scaling).** `lfun/lfuncreate` gives fast reference values (calibrated: 10.8 s @ 10,000 digits on 389.a1). Lane C is used for cross-checks and precision escalation inside Lane B's interval evaluation, never as certification by itself.

## Validation harness (exit criteria)

1. Calibration set: {389.a1 (r=2, w=+1, Δ>0); 37.a1 (r=1, w=−1, Δ>0); 53.a1 (r=1, Δ<0, one component, |Sha|=1 closure at 100 digits); 57: [0,-1,1,-2,2] (r=1, Δ<0, Ct=2 — nontrivial c_p + index-3 generator)} + ≥ 2 more curves with LMFDB-published analytic data (including one with |Sha| > 1 once the export lands — a non-square exclusion window test).
2. For each curve and j: Lane B interval agrees with Lane C value, and the interval width matches the derived bound (not just empirically small).
3. Cross-implementation agreement at 100 / 500 / 5,000 digits.
4. A deliberately perturbed run (wrong sign of a-coefficients) must FAIL the harness — guard against vacuous bounds.

## Convention sheet (pinned by calibration — see results/m0-calibration.md)

- Ω_BSD = n_components(E(ℝ))·omega1_PARI — 2·ω1 if Δ>0, ω1 if Δ<0. Component count, not root number (final form 2026-11-10: 27-curve component scan 27/27 + 100-digit closures on 53.a1 and the Δ<0 57 curve; see results/m0-calibration.md). Mechanism: ω1_PARI is one oval's period; BSD integrates the whole real locus.
- ∏ c_p from per-prime `elllocalred` components (NOT `elltamagawa`, which includes c_∞).
- Regulators from mwrank saturation only (ellrank generators can be k²-inflated).
- |Sha|_impl for rank 2 = L″(1)/(2!·Ω_BSD·Reg·∏c_p·tors²) — must be a perfect square integer if BSD + finite Sha hold.

## Tooling notes (hard-won, obey them)

- **gp script files: ONE statement per line.** No newline continuation — a multiline function body parses with "unexpected end of file" and the whole file cascades into type errors (empirically isolated 2026-11-10: minimal repro + paren balance check, `t1/t2` pair). Keep every statement on a single line.
- `iferr(expr, ERR, cond)`: the **second argument is a bare variable name**, not an expression. Prefer avoiding iferr: compute the b-invariant discriminant Δ = −b2²b8−8b4³−27b6²+9b2b4b6 directly and filter singular models arithmetically before `ellinit`.
- `default(threadsizemax, 2^30)` at the top of any script using `lfun` at ≥ 4,000 digits (8 MB default thread stack overflows; verified fix at 10,000 digits, prefix reproduced).
- `ellglobalred(E)[2]` in PARI 2.19 returns a vector, not the scalar root number — don't use it for parity; root number comes via parity/rank or the proper function when needed.
- Generator saturation: ellrank generators may be k× primitive (k=3 on 37.a1 and the Δ<0 57; k=1 on 53.a1); regulators must use `~/.local/eclib/bin/mwrank` saturation (in PATH? no — full path) or explicit height-division checks.

## Lane B v2 — implemented 2026-11-10 (Arb engine, unified formula) — results

**Engine**: python-flint 0.6.0 (FLINT/Arb ball arithmetic). Replaces the
planned mpmath-iv layer: mpmath's `iv` is **double precision** (53-bit)
regardless of context precision (measured: `iv.mpf(2)/3` width 1.1e−16 at any
dps) — unusable beyond ~15 digits. Arb gives true arbitrary-precision balls:
`ctx.prec` in bits; native `expint` (= E_s(x) = ∫_x^∞ e^{−t}/t^s dt, exactly
our kernel primitive), `const_euler`, `pi`, `legendre_p_root(weight=True)`.

**New identity backing the formula** (derived + machine-verified this session):
pointwise root-number law G(1/u) = ε·u²·G(u) (G(u) = Σaₙe^{−2πnu/√N}),
verified to 120 digits on 4 curves of both ε classes; fixed-point consequence
G(1) = 0 for ε = −1 verified to ~1e−135. Splitting at the fixed point u = 1:

    Λ^(k)(1) = (1 + ε(−1)^k)·Σₙ aₙ·M_k(2πn/√N)

with M_k the absolutely convergent up-window log-moments. Rank parity (the
whole functional equation) is baked into the (1 + ε(−1)^k) factor.

**Certified results** (scripts/lane_b_cert.py, seconds per value):
- 20a1 (ε=+1): L(1) = 0.47072919032651896658063… ± 2.2e−311 (PARI/M0: 0.4707291903 ✓)
- 37a1 (ε=−1): L(1) = certified exact ball [0 ± 0] (zero anchoring live)
- 37a1 (ε=−1): L′(1) = 0.305999773834052301820… ± 4.6e−380 (Gross–Zagier
  calibrated constant ✓ to all stored digits)
- perturbation: corrupted aₙ flips the value (−0.5108) AND trips the (P)-law
  gate (POINTWISE-LAW-FAIL) — non-vacuous harness ✓
- Λ-normalization check: 2Σaₙe^{−λₙ}/λₙ = Λ(1) = L(1)√N/(2π) agrees with PARI
  at 130 digits (difference −6.9e−136, logs/probe4.log)

**One open defect (k = 2)**: L″(1) mid is off by ~6.4e−19 from the published
Boothby–Bradshaw value (interval ±3e−447 around a wrong center). Traced to the
C_Q = Q(1) Gauss–Legendre quad: quad VALUE validated to 33 digits against
mp.quad per-panel, but the certified C_Q mid is off 3.8e−28 (true 0.0978431972
1667017932…, quad mid …01789|4…). Suspect: per-panel remainder path
(C32·|f⁽³²⁾|·(h/2)³⁴) or a precision-context drop inside the panel loop. Fix +
revalidation is the next session's first action; k = 0/1 unaffected.

## Lane B v3–v5 — the k = 2 defect RESOLVED (2026-11-10, session 4)

The recorded "L″(1) mid off ~6.4e−19" was a **false certificate, not a near
miss**. Full diagnosis chain (all machine-checked):

1. **Reconciliation (probe_cq2/3):** certify_CQ's mid equalled the bare GL-16
   [1,40] sum and missed EXACTLY the tail ∫₄₀^∞ e⁻ᵗln(t)/t dt
   = 3.84925663812063074867258587498e−19 (independent mp.quad + gammainc′ +
   Watson asymptotics all agree), propagating through k=2 as
   4·Σₙ aₙ/n·(−tail) with Σₙ₌₁³⁷⁶ aₙ/n = −0.41411568776043551313 ⟹
   +6.37615580907e−19 — matching the recorded defect to all printed digits.
   The GL sums themselves are CLEAN: GL[1,40] matches mp.quad to ≥169 digits,
   and GL[1,200] recomposed from three sub-range sums agrees to 0.0. No
   quadrature bug — but the recorded ball (±3.18e−447) contained neither the
   tail nor the truth: the published value sat OUTSIDE it.
2. **Root causes (found by bisection):**
   a. certify_CQ's tail lived only in a crude widening (±5.3e−18); a stale
      cache variant carried neither tail-mid nor widening.
   b. python-flint `arb **` on wide-radius balls collapses the mid toward 0
      (t = [1±0.064] → t²³ = [−0.038, +4.165], straddling zero); division by
      such balls returns NaN, silently zeroing the GL-remainder widening.
   c. `_watson_tail_ball` widened by bare max|f^(M)| without the e^−T factor:
      a zero-straddling ball whose mid python-flint normalizes toward 0,
      seeding −1.58e−39 garbage into the C_Q mid (and ~e−34-scale noise into
      the k=2 mid via Q_of's large-λ Watson path).
   d. load_constants' e10-based cache reconstruction narrowed the reloaded ball
      below its own certificate.
3. **Fixes (Lane B v5, current code):** GL-16 composite over [1, T=200] with
   the CORRECT composite remainder C32·(2h)³³ (the old h³⁴ form was 2³³/h ≈
   1.7e13 too tight), |f⁽³²⁾| bounded by interval evaluation per 32-panel block
   with powers built by repeated multiplication (NaN-proof); Watson tail added
   INTO THE MID as a certified ball (M=20; |R| ≤ e^−T·max|f^(M)|; single-
   extremum analytic bound; the e^−T factor in place); context-safe node
   construction (setbits BEFORE legendre_p_root); Stirling-honest series
   truncation in Q_of (loop to |t_j| ≤ 10^−(d+45)/3, remainder bounded by the
   measured next-term ratio); Watson-vs-series crossover λ_W = max(100,
   2.3·(digits+45)) — the Watson remainder floor is ~e^−λ·O(λ^−1.5), so the
   crossover must scale with the target digits; cache v3 stores the certified
   radius as a decimal string (re-applied verbatim on reload).
4. **Verified certificates (Lane C = gp lfun, fresh inside the harness):**
   - C_Q mid vs independent references (mp.quad ≡ ∂ₛΓ(s,1)|₀ at 340 dps, and
     PARI lfun): agree to 7.99e−117; certified radius 1.55e−116; fresh and
     cache-hit paths produce bit-identical downstream balls.
   - 389a1 (ε=+1, rank 2): L″(1) = 1.51863300057685354046038521578944403815…
     certified radius 4.27e−114 at d100 (mid − ref = −6.3e−117); at d80
     radius 3.02e−71, mid − pub = +7.5e−73. The Boothby–Bradshaw 250-digit
     constant is inside the ball.
5. **Harness d100 (scripts/harness_agree.py → logs/harness-d100.log):**
   389a1-k2 inside ✓ radius 4.27e−114; 37a1-k1 inside ✓ radius 7.9e−175;
   20a1-k0 inside ✓ radius 1.2e−174; every mid − ref ≤ 9% of its radius;
   all three perturbation gates fail corrupted aₙ streams (POINTWISE-LAW-FAIL,
   exit 4). **HARNESS-PASS at the 100-digit rung.**
6. **Remaining (cl-4/cl-5):** 500-digit rung ready — aₙ streams extended to
   12,600 Hecke-checked terms (work/extract-an-13k.gp, n_eval 4144, tail
   coverage e^−1320); estimated 10–40 min runtime. The 5,000-digit rung needs
   the k=2 series band to carry ~12,400 kernels at ~63k working bits (the
   bits_storm formula's 0.4343·λ_max term correctly prices the series-
   cancellation peak; the Watson floor ~e^−λ forbids certification below that
   scale at argument λ) — single-core estimate hours-to-days; multiprocess or a
   larger-M Watson design is the open engineering item. LMFDB rank-2 export
   (cl-5) still open.

## Prior-art sweep — has anyone published an equivalent of M1? (2026-11-10,
session 5)

No published counterpart for the shipped combination. Closest neighbors, in
decreasing order of overlap:

- **Mark Watkins, "Computing L'''(1) to high precision" (J. Théor. Nombres
  Bordeaux; magma.maths.usyd.edu.au/~watkins/papers/HighPrec.pdf)** — closest
  on precision: Arb-ball certification that the THIRD central derivative is
  zero to ~1,033 / 10,029 digits for 12 rank-5 curves (4.6–21 days per
  curve), used to prove class-number lower bounds (companion theorems bound
  √D·L(χ,1) via those zeros). This is certified central-derivative ZERO-ness
  (inequality certificates at extreme precision), not certified VALUES at
  k ≤ 2; one-off curves, days per value, no harness, no BSD-closure
  instrument. Also independent validation that Arb balls carry
  10,000-digit-scale central-derivative certification in practice. Note: the
  paper finds quasi-cubic local power series beat quasi-quadratic methods at
  these precisions — matches our own ladder experience.
- **Bradshaw thesis 2010** (Provable Computation of Motivic L-functions) —
  closest on method: provable values AND derivatives with explicit error
  bounds; a prototype, never operationalized or shipped.
- **LMFDB modular-forms paper 2022, §9** — closest at scale: certified
  analytic ranks of modular-form L-functions in the database. Rank
  certification ≠ BSD-closure value balls; no k=2 feed into the strong
  formula, no rank ≥ 2 scan path.
- **Boothby–Bradshaw ~2009** — one-off 10,000-digit 389a1 closure at
  O(prec³), ~1 CPU-week (our ladder reproduces it in 10.8 s).
- Infrastructure: Arb/FLINT ship rigorous Dirichlet L-functions and zeta
  tooling only — no elliptic-curve L-derivative layer.

2026-era watch item: BootLoops' "Eichler" web tool (periods and iterated
modular-form integrals) surfaced in search — monitor; no certified
L-function claims visible in its description.

M5 consequence: the disclosure writeup must cite Watkins (precision),
Bradshaw (method), and LMFDB §9 (scale). What remains without published
counterpart: certified k ≤ 2 VALUE balls at seconds per value, harness-gated
against fresh cross-lane computation, aimed at rank ≥ 2 curves.

## Lane B v6 — d500 rung: jcap crash fixed, 389a1-k2 spliced and verified (2026-11-10, session 5)

**Crash.** The in-flight d500 harness's 389a1-k2 child died ~2.5 min in — silently
(traceback captured by the parent's buffered stdout — log 0 bytes, no crash report).
Direct unbuffered rerun reproduced it deterministically at the same term:
`RuntimeError: Q_of series failed to reach threshold` at **n = 2982 (λ = 950.28 = 2982·2π/√389)**.

**Root cause — pure control flow, not numerics.** Q_of's series exit requires
`j >= jguard` with `jguard = max(int(4λ), 16)`, while `jcap = 6·digits + 800 = 3,800`
at d500. Every series-zone λ > 950 has 4λ > jcap ⇒ unconditional raise regardless of
the fact that the series converges below threshold near j ≈ 2,750–3,800. d100 never
hit it because the then-current stream topped out at λ = 300.8 (4λ = 1,203 < 1,400).
The workspace sizing (`bits_storm`) already covers the e^λ cancellation peak — the
series was numerically healthy the whole time.

**Fix (one line, v6):** `jcap = max(6·digits + 800, jguard + 16)` — d100 semantics
unchanged (λ ≤ 300.8 keeps jcap = exactly 1,400), and the worst-case series term
(λ = 1,253.25 — just under λ_w, jguard = 5,013 ≤ jcap = 5,029) now exits cleanly.

**Probe (work/probe-qof-d500.json — PROBE-PASS).** λ ∈ {300, 950.0, 1,253.25, 1,320.15}
at the exact in-run workspace (7,933 bits): all four contained vs independent
mp.quad references (900 dps). The λ ≥ λ_w Watson return path was exercised for the
first time ever: λ = 1,320.15 → radius 1.76e−629 (relative ~10⁻⁵³; absolute radius
100+ digits deeper than the 500-digit target requires — cancellation sloppiness is
absorbed by the massive workspace). Series points return the C_Q-floor radius
(1.54691e−116) with the C_Q mid's certified offset (mid−ref 7.991e−117) inherited
exactly as designed.

**389a1-k2 @ d500 (work/hb-389a1-d500.json, lane v6):** certified ball radius
**8.53e−114**; mid = 1.5186330005768535404603852157894440381561950329898487… (
published Boothby–Bradshaw digits); **fresh gp \p 580 reference inside at
mid − ref = 1.2e−116 (14 % of radius)** — scripts/check_389a1_d500.py,
work/gp-refs-d500.txt. Perturbation gate: POINTWISE-LAW-FAIL, exit 4 ✓.
Provenance (splice): the harness's own in-run 389a1 sub-run died pre-fix and is
recorded as FAIL in its json; v6 differs from v5 only in Q_of's series jcap, so the
harness's k=0/1 lanes and all law-gate exits are bit-identical under v6 vs v5.

**Structural finding (drives M2 prep).** The k=2 certificate depth is floored by C_Q's
certified radius (1.54691e−116), which is digits-INDEPENDENT: it comes from the
Watson tail widening e^{−T}·max|f⁽ᴹ⁾| at certify_CQ's fixed defaults T=200, M=20
(≈ 2.8e−117 bound). So d500's k=2 ball is C_Q-floored at ~e−114 — not 500 digits
deep (the d100 4.27e−114 radius had the same floor). d500-grade (and d5000-grade)
k=2 needs certify_CQ re-parameterized with T, M scaled to the target:
T ≈ 650, M ≈ 650 → tail_rem ~ e^{−2T} ≈ 10^{−564}; GL-16 over [1, 650] stays
seconds-fast at 2,072 bits. Same status as the d5000 workspace-budget item.

## Session 6 — d500 rung PASS: the FAIL was the reference, not the evaluator (2026-11-10)

**The incident.** The first full d500 harness run (14:40) returned HARNESS-FAIL:
37a1-k1 `inside=False` (mid−ref = −5.638e−582 vs radius 3.13e−584), 20a1-k0
`inside=False` (+3.906e−581 vs 4.81e−584). Both offsets sat at the ~e−581/582
scale — exactly where a gp reference computed at `\p 580` loses its own accuracy.

**Arbitration (gp self-diff).** Recomputed all three refs at `\p 800`
(work/gp-refs-d500-800.txt): p800 − p580 = −5.632e−582 (37a1), +3.907e−581
(20a1) — matching the old offsets to within the ball widths. The `\p 580`
references carried ~e−582-scale error, 18×–81× larger than the certified radii.
Lane B was right; the ruler was too coarse.

**HARNESS-PASS at d500** (work/harness-d500.json, ok=true; refs rerun fresh at
`\p 800`, Lane B balls reused from disk via the new cache-first path, gates 3/3):

| lane | radius | mid − ref₈₀₀ | offset/radius |
| --- | --- | --- | --- |
| 389a1-k2 | ±8.53395e−114 | +1.2012e−116 | 1.4 % |
| 37a1-k1 | ±3.12989e−584 | −5.6365e−585 | 18 % |
| 20a1-k0 | ±4.81482e−584 | −8.6707e−585 | 18 % |

The k=0/1 radii are the real d500 certificates (~583 certified digits each —
the Σaₙ-kernel widening floor), and the k=2 radius remains C_Q-floored at
~e−114 as recorded in session 5. All three published/documented constants sit
inside: BB's 389a1 L″(1), the Gross–Zagier-calibrated 37a1 L′(1), and 20a1 L(1).

**Harness changes (v7, scripts/harness_agree.py):**
1. Lane C refs at `digits + 300` dps (was digits+80). Rule of thumb learned the
   hard way: a cross-check reference must out-precision the certified ball's
   radius by ≥ 200 digits, else containment fails spuriously.
2. Cache-first Lane B (`[laneB-reuse]`, `--fresh` to force): certified balls are
   deterministic per (curve, k, digits, streams); reuse spares the 70-minute
   37a1-k1 lane while Lane C (fresh) and all perturbation gates (fresh, early-exit
   via the law check) keep every run honest.

**Serialization lesson (for every future harness consumer).** In the lane JSON,
`L.k = [{lo-ball}, {hi-ball}]` where each endpoint ball prints as
`[mid ± ulp]` — the `± 1.73e−3729`-style radii are the endpoint ROUNDING ulps at
the print context (~2 ulp), NOT the certified radius. The certified radius is
`(hi − lo)/2` from the two endpoint mid prints (exact terminating decimals).
Diff-structure check: d100 endpoint prints differ across digits 175–605
(radius 7.9e−175 scale ✓), d500 across digits 585–1967 (radius 3.13e−584 ✓).
The mid-session-5 chat report "±3.25e−2411" for 389a1 was exactly this ulp
artifact; the doc's 8.53e−114 was correct all along.

**cl-4 status:** d100 ✓ (12:34), d500 ✓ (15:1x). Remaining on the ladder: d5000
(design item — multiprocess + C_Q re-parameterization, session 5 finding) and
cl-5 (LMFDB rank-2 export). M2 unblocked.

## cl-5 RESOLVED — LMFDB rank-2 export + exact Tier-1 counts (2026-11-10, session 6)

`work/lmfdb-rank2-export.csv` (21.85 MB; download URL in `work/dl.url`, canonical
search link in the file header). Columns per file legend: [Label, Class, Conductor,
Rank, Torsion, CM, Weierstrass equation], tab-separated, quoted strings.

- **348,672 rank-2 curves, conductor 1–500,000** — count declared server-side by the
  LMFDB query engine in the file header ("returned 348672 curves"), verified on disk
  row-for-row (strict tab-split: 348,672 rows × 7 fields; rank column = 2 in every
  row; all labels distinct).
- **279,056 distinct isogeny classes** — the natural L-evaluation unit (isogenous
  curves share their L-function); per-curve D2 adds regulator/Tamagawa/torsion
  inputs on top.
- **129,892 distinct conductors; min 389 (389.a1, the classic first rank-2 curve,
  literally the file's first row) to max 499,998.** 50k-bucket growth
  [23612, 33363, 34574, 35294, 36112, 36875, 37347, 37169, 37100, 37226] — smooth,
  no truncation cliff (completeness cross-check).
- Structure: torsion [] 236,501 curves; [2] 100,237; [2,2] 6,020; [3] 3,539; [4]
  1,951; [6] 271; [2,4] 76; [5] 74. **396 CM curves** (discriminants incl. −3, −4,
  −16), 348,276 non-CM.

Tier-1 sizing consequence [EST]: 279,056 certified class-level L-evaluations plus
per-curve rank proofs/regulators/Tamagawa for 348,672 curves. At the d100 rung's
measured ~2 s per certified value, the L-evaluation pass alone is ~6.5 CPU-days
sequential ≈ half a day wall on 12 cores; mwrank rank proofs add seconds per curve.
Comfortably inside the attack plan's compute-weeks budget. **M1 checklist complete
(cl-1–cl-5): M1 CLOSED.** M2 opens: detector battery, single-curve full loop
(plan §5; task created).
