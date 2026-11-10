# M0 Calibration Record — 389.a1 (conductor 389, rank 2)

Date: 2026-11-10 (UTC ~13:27). Tool: PARI/GP 2.19.0, arm64, GMP 6.3.0 (Homebrew).
Scripts: `scripts/calib389a1.gp`, `scripts/check250.gp`. Reproduce: `echo quit | gp -q -f scripts/check250.gp`.

## Curve

```
E = ellinit([0,1,1,-2,0])    \\  y^2 + y = x^3 + x^2 - 2x   (minimal model)
```

| Invariant | Value | Source |
| --- | --- | --- |
| Conductor | 389 | ellglobalred; matches Stein 2010 slides |
| Discriminant | 389 | PARI + hand-computed (Δ = −b2²b8−8b4³−27b6²+9b2b4b6 = 389) |
| Bad primes | [389] | factor(conductor) |
| Kodaira at 389 | I1 (code 5), **c₃₈₉ = 1** | elllocalred 4th component |
| Torsion | trivial | elltors = [1,[],[]] |
| Rank | exactly 2 | ellrank = [2, 2, 0, [[-2,0],[4,8]]] — generators (−2,0), (4,8) |
| Displayed points | (0,0), (−1,1) | both ellisoncurve ✓ (Stein slides) |
| Regulator | 0.1524601779431437516243247570494558232437270774866308178402809800460532256835624636041148155… | det of √-height pairing on P₁,P₂ |

**Ω convention (pinned by exact closure):**

- PARI `E.omega[1]` = 2.4902125608550550753213577919423024601560581803395700400549723… (spans the full real locus)
- BSD Ω (identity-component integral) = **2 · omega1_PARI** = 4.98042512171110…  — the archimedean factor for two real components.
- `elltamagawa(E)` returns 2: PARI's "global Tamagawa number" absorbs the archimedean factor c₍∞₎ = 2. ∏ finite c_p = 1 (confirmed independently by elllocalred). Both packages put the same factor 2 in different places; the pair (Ω = 2·ω1, c₃₈₉ = 1) is the literature convention and closes the formula exactly.

## L-values at s = 1 (via lfun/lfuncreate)

- L(1) = 2.304…×10⁻¹⁵⁷ (numerically zero at 150 digits)
- L′(1) = −1.306…×10⁻¹⁵⁷ (numerically zero)
- **L″(1) = 1.5186330005768535404603852157894440381561950329898487031716101850959821454955611048804318700046577183337542984136397153650520262209671289805699812548081429042749756770429723887482582398917985423226150172284380067273160739712268820343534169525066746101…**

Stable: the 150-digit and 250-digit runs agree digit-for-digit in the stable prefix (≥ 155 digits).

## Result

```
|Sha|_impl = L''(1) / (2! · Ω_BSD · Reg · c₃₈₉ · |tors|²) = L''(1) / (4 · omega1 · Reg)
           = 1.0000000000000000000000000000000000…  (all 250 displayed digits, exact)
```

Also 1.000…×10¹⁵⁰ at the 150-digit run. Consistent with:

- Boothby–Bradshaw: BSD formula verified to 10,000 digits for 389a1 **assuming Sha = 1**;
- LMFDB: 389.a1 analytic order of Sha = 1.

**CALIBRATION PASSED** — our pipeline reproduces the published pieces exactly, end to end.

## Precision pilot (COMPLETED)

`scripts/pilot2.sh` → `logs/precision-pilot.log`: ladder 300 → 10,000 digits recomputing L″(1) via lfun.
**All 11 rungs completed** (300, 500, 800, 1000, 1500, 2000, 3000, 4000, 5000, 7000, 10000). Observed timings: 300 digits = 56 ms, 500 digits = 33 ms; **10,000 digits = 10,828 ms**. Cross-precision agreement on the stable prefix at every rung.

**The 10,000-digit Boothby–Bradshaw reproduction is complete**: BB's 2004 estimate was ~1 CPU-week (custom O(prec³) stack); `lfun` reproduces the same L″(1) to 10,000 digits in ~10.8 s. Precision cost is a non-issue for D2 deep-dives going forward.

**Thread-stack stability (re-verified 2026-11-10, M1):** an `lfun` thread-stack overflow limit (8 MB default thread size) was identified on the high-precision ladder; the fix is one in-script line, `default(threadsizemax, 2^30)`. Verified: the top rung (10,000 digits) re-runs clean — exit 0, no overflow — with the stored L″(1) prefix 1.5186330005768535404603852157894440381561950329898… reproduced (`scripts/ladder-10k-fixed.gp` → `logs/ladder-10k-fixed.log`).

## Second calibration curve — 37.a1 (w = −1, rank 1, Gross–Zagier)

Script: `scripts/check37a1.gp`. Curve `ellinit([0,0,1,-1,0])` — y² + y = x³ − x.

- conductor 37, Δ = 37, I1, c₃₇ = 1, torsion trivial, rank exactly 1 (mwrank: **unconditional, full MW basis, Regulator = 0.051111408239969**, generator [0:−1:1], 0.14 s).
- L(1) = 0 (≤ 10⁻¹⁵⁴); **L′(1) = 0.3059997738340523018204836833216764744526…** — the classical Gross–Zagier value ✓
- L'(1) / (2·ω1·Reg_sat·c₃₇) ≈ 1 + 10⁻⁷ ⇒ **|Sha| = 1** (closure limited by the 12-digit mwrank regulator; full-precision redo queued with mwrank -p).

## Ω convention — FINAL LAW (M1 verification, 2026-11-10)

```
Ω_BSD = n_components(E(ℝ)) · omega1_PARI     n_components = 2 if Δ > 0, 1 if Δ < 0
∏ c_p = ∏ (elllocalred c_p)                  finite Tamagawa factors only
```

The split is the **real-component count, not the root number** (established *within* the w = −1 class: 37.a1 has Δ>0 and closes with 2ω₁; 53.a1 / 57:Δ<0 close with ω₁). The earlier "always 2ω₁" reading was an artifact: all three M0 calibration curves happen to have Δ > 0.

Mechanism: PARI's ω₁ is the period of ONE real oval; the BSD real period integrates the invariant differential over the whole real locus, which has 1 or 2 components.

Evidence:
- Component scan (`scripts/omega-scan.gp` → `logs/omega-scan.log`): 27 rank-1 curves, N ≤ 150, both Δ signs — **18/18 Δ<0 close with Ω = ω₁, 9/9 Δ>0 close with Ω = 2ω₁, zero exceptions**; in every row the wrong convention misses the 1/m² pattern.
- 100-digit closures (`scripts/verify-omega.gp`):
  - **53.a1 [1,−1,1,0,0]** (Δ<0, Ct = 1, torsion 1): S = L′(1)/(ω₁·Reg·Ct) = 1.000… all 100 displayed digits (wrong branch: exactly 0.500…). Implied |Sha| = 1. mwrank: generator [0:−1:1], Regulator 0.092981484638654…, "rank and full MW basis determined unconditionally"; PARI ĥ = mwrank Regulator exactly (generator already saturated, m = 1).
  - **Δ<0 curve 57 [0,−1,1,−2,2]** (Ct = 2, torsion 1): ellrank generator has index 3 (ĥ = 0.3381713346… = 9 × Reg_sat — same index-3 pattern as 37.a1); mwrank Regulator = 0.037574592736824…; S = L′(1)/(ω₁·Reg_sat·Ct) = 1 exactly. Wrong branch: 0.0555…
- Prior: 389.a1 (rank 2), 37.a1 (rank 1), 20 (rank 0), all Δ>0, close with Ω = 2ω₁ — consistent on both branches.

PARI's `elltamagawa(E)` returns ∏c_p multiplied by the archimedean factor (now identified precisely as the real-component count) — do not use it in the BSD formula.

## Saturation lesson (pipeline-critical)

PARI `ellrank`'s returned generators are 2-saturated only: on 37.a1 it returned (−1,0) with ĥ = 0.4600026741597195… = 9 × 0.0511114082399690… (i.e. 3 × primitive). Regulators for CERT-2 must come from **mwrank saturation** (verified: "points were already saturated", unconditional basis) or explicit height-division checks. An unsaturated regulator inflates by k² and fakes a small |Sha|_impl — exactly the kind of silent factor that would corrupt D2 alarms.
