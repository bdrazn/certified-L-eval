# certified-L-eval

Certified interval evaluation of **L(1), L′(1), L″(1)** for elliptic curves E/ℚ — machine-checked balls with derived error budgets, seconds per value at the 100-digit rung, harness-gated against fresh PARI/GP computation at the 100- and 500-digit rungs.

Part of an open campaign to attack the Birch–Swinnerton-Dyer conjecture with machine-checkable certificates: [docs/01-survey-bsd-state-of-art.md](docs/01-survey-bsd-state-of-art.md) (state of the art) · [docs/02-attack-plan.md](docs/02-attack-plan.md) (plan) · [docs/m1-certified-evaluator.md](docs/m1-certified-evaluator.md) (this tool's design + full validation log) · [docs/m2-detector-battery.md](docs/m2-detector-battery.md) (detector battery: full-loop audits + validation batch).

## What "certified" means here

Every output is a ball `[lo, hi]` whose width comes from a **derived** budget — Hasse tail majorant (`|aₙ| ≤ 2√n`), quadrature remainders, rounding via Arb ball arithmetic — not from empirical agreement. The engine is [python-flint](https://github.com/flintlib/python-flint) (FLINT/Arb, true arbitrary-precision balls; mpmath's `iv` was measured double-precision-only and is unusable past ~15 digits).

The core identity (derived and machine-verified, re-checked on every run — not assumed): for G(u) = Σ aₙ e^(−2πnu/√N) on a strong Weil curve, the pointwise involution law G(1/u) = ε·u²·G(u) (ε = root number) splits the Mellin transform at its fixed point u = 1 into

```
Λ^(k)(1) = (1 + ε(−1)^k) · Σₙ aₙ · M_k(λₙ),     λₙ = 2πn/√N
```

with absolutely convergent kernels M₀ = e^(−λ)/λ, M₁ = E₁(λ)/λ (Arb `expint`), M₂ via a once-certified constant C_Q (GL-16 composite quadrature + certified Watson tail). Rank parity is baked into the prefactor: **ε = −1, k = 0 returns the exact ball [0, 0]** (zero anchoring), and 37a1's L′(1) > 0 is proven rigorously against a closed-form-calibrated Gross–Zagier constant.

Before any evaluation the stream must pass selftest gates: Hecke multiplicativity checks (exit 3 on violation) and the pointwise involution law (exit 4 on violation). A deliberately sign-corrupted aₙ stream must FAIL — the gates are tested against corruption, not assumed vacuous.

## Validation (harness-gated)

`scripts/harness_agree.py` re-derives Lane C = fresh PARI/GP `lfun` values inside each run and demands containment: the reference must lie **inside** the certified ball. Artifacts: [certificates/](certificates/).

| lane | d100 ball | d500 ball |
| --- | --- | --- |
| 389a1 · L″(1) · r=2, ε=+1 | ±4.27e−114 · mid−ref −6.3e−117 | ±8.53e−114 · mid−ref +1.2e−116 |
| 37a1 · L′(1) · r=1, ε=−1 | ±7.88e−175 · mid−ref +7.2e−176 | ±3.13e−584 · mid−ref −5.64e−585 |
| 20a1 · L(1) · r=0, ε=+1 | ±1.21e−174 · mid−ref +1.1e−175 | ±4.81e−584 · mid−ref −8.67e−585 |

- Every fresh reference lies inside its ball; every mid−ref offset ≤ 9% (d100) / 18% (d500) of its radius.
- All three perturbation gates caught corrupted streams at both rungs (exit 4, POINTWISE-LAW-FAIL). **HARNESS-PASS at d100 and d500** (`certificates/harness-d100.json`, `certificates/harness-d500.json`).
- The certified values agree with the literature constants digit-for-digit: 389a1 L″(1) = 1.518633000576853540460385…, 37a1 L′(1) = 0.305999773834052301820…, 20a1 L(1) = 0.47072919032651896658063…
- Speed: ~2 s per certified value at d100. The d500 rung ran ~1–2 h per lane on a 12-core M-series Mac. The 300→10,000-digit precision ladder on 389a1 reproduces the Boothby–Bradshaw 10,000-digit closure end to end, 10,000 digits in **10.8 s** (their 2004 estimate for the same: ~1 CPU-week, O(prec³) method).

## Field-tested pitfalls

Why the harness exists — four ways the tool's predecessor certified a wrong number, all caught and documented (full chain in [docs/m1-certified-evaluator.md](docs/m1-certified-evaluator.md)):

1. **A wrong ball that looks certified.** The first k=2 certificate missed a quadrature tail (∫₄₀^∞ e^(−t)ln(t)/t dt = 3.84925663812063074867258587498e−19) and printed a ball ±3.18e−447 around a center off by +6.37615580907e−19 — the published truth sat *outside* the "certified" interval. Found by reconciling mid−ref at extreme precision; the quadrature itself was clean to ≥169 digits.
2. **python-flint `**` on wide-radius balls collapses the mid toward 0**; downstream divisions silently NaN-poison the widenings.
3. **Mid-normalization on zero-straddling balls** seeds garbage mids into widenings.
4. **References must out-precision the ball.** At 500 digits the ball radii (±e−584) are *tighter* than PARI's own self-error at `\p 580` (~e−582) — containment then fails spuriously. Measured gp self-diff p580→p800: 5.6e−582. Rule pinned in the harness: refs at `\p digits+300`.

## Detector battery (v0.2.0)

`scripts/detector.py` turns the evaluator into a per-curve auditor. One command runs the full loop — Hecke-checked aₙ stream (v7 conductor-scaled sizing) → mwrank unconditional rank + regulator → PARI root number / Tamagawa / torsion / ellrank bounds → certified L^(k) balls (v7) → fresh out-precision `lfun` references → containment gate — and stamps every number's provenance (command, exit code, elapsed, sha256 of every file). Then three detectors demand agreement:

- **D1** certified non-vanishing ⇒ ord ≤ k (k matched to the root number)
- **D3** parity: proven rank parity vs root-number parity
- **D2** implied-|Sha| square test (rank 2): certified interval for `|Sha| = L″(1)·T² / (2!·Ω·Reg·∏c_p)` must contain a perfect-square integer

Seeded validation batch (reproducible: `--sample 10 --seed 20261110`; 389a1 pilot included): **11/11 CLEAN** — all containments pass, mwrank ranks agree with PARI bounds, regulator cross-checks at machine precision, and every implied-|Sha| window contains exactly [1]. The batch itself caught and fixed two real defects (v6's small-N-only selftest sizing; a D2 torsion-factor placement bug that fired the battery's first two alarms on torsion-2 curves — post-fix they read |Sha| = 1.000… to 38 digits). Full story: [docs/m2-detector-battery.md](docs/m2-detector-battery.md); artifacts in [certificates/m2/](certificates/m2/).

Runtime: 2 s (N = 389) to ~55 s (N ≈ 500k) per full audit at d100.

## Quickstart

```bash
pip3 install python-flint mpmath          # + PARI/GP ≥ 2.19 for Lane C references

# single certified value
python3 scripts/lane_b_cert.py --N 389 --an work/an-389a1.txt --eps 1 --k 2 --digits 100

# the full gated harness, 100-digit rung
python3 scripts/harness_agree.py --digits 100        # add --fresh to recompute all lanes
python3 scripts/harness_agree.py --digits 500        # needs the 12,600-term streams (included)

# the M2 detector battery: full per-curve audit (rank, parity, certified L, |Sha|)
python3 scripts/detector.py --label 389.a1 --digits 100
python3 scripts/detector.py --sample 10 --seed 20261110 --digits 100
# needs also: eclib's mwrank at ~/.local/eclib/bin/mwrank; the LMFDB rank-2 export at
# work/lmfdb-rank2-export.csv (details: docs/m2-detector-battery.md)
```

`work/` ships the validated aₙ streams (12,600 Hecke-checked terms per curve, 389a1/37a1/20a1) and the certified constants cache. Harness outputs land in `work/hb-*.json` (gitignored); the shipped validated runs are archived under `certificates/`.

## Status — honest

- **Validated:** everything in the tables above, at two rungs, against two independent implementations, with corruption-tested gates.
- **Not yet:** this tool reproduces constants the literature already holds — it is an instrument, not a result. The claim that "certified k ≤ 2 *value* balls at seconds per value, harness-gated, aimed at rank ≥ 2 curves" has no published counterpart *as far as our own sweep went* (see the prior-art section of [docs/m1-certified-evaluator.md](docs/m1-certified-evaluator.md)) — that sweep is ours, not referee-grade.
- **Roadmap:** audit-trail pipeline + detector battery **shipped (v0.2.0)**; Tier-0 bit-agreement bridge **passed 100/100, zero alarms** (seeded random rank-2 sample at d100, mean 39.4 s/curve — [certificate](certificates/m2/bridge-tier0-seed20261110-d100.json)). Tier-1 verification tables over the LMFDB rank-2 family, conductors < 500,000 (exact counts, verified: **348,672 curves / 279,056 isogeny classes / 129,892 conductors**). Tier-0-corrected sequential cost: **≈13.7M CPU-seconds ≈ 159 CPU-days at d100** (bridge mean 39.41 s/curve, range 8.7–173.4 s; the earlier "2–55 s ⇒ a few CPU-days" estimate under-measured mid-family curves) — the scan was launched in parallel (8 workers over a conductor-sorted queue) and **paused on 2026-11-12 at 1,062 curves — all clean, zero alarms; frontier conductor ≈ 6,354; fully resumable** (parked snapshot: [certificates/m3-tier1/parked-snapshot-20261112.json](certificates/m3-tier1/parked-snapshot-20261112.json); resume with `python3 scripts/tier1_scan.py --supervise --workers 8`). No preprint before those tables exist.

## Prior art & credits

- R. Bradshaw, *Provable Computation of Motivic L-functions* (PhD thesis, 2010) — provable values + derivatives with explicit error bounds; a prototype, never operationalized.
- M. Watkins, *Computing L'''(1) to high precision* (J. Théor. Nombres Bordeaux) — Arb-ball certification of central-derivative zero-ness at 1,033–10,029 digits for 12 rank-5 curves.
- Boothby–Bradshaw (~2009) — the 10,000-digit 389a1 closure our ladder reproduces.
- LMFDB modular-forms paper (2022, §9) — certified analytic ranks at database scale.
- Infrastructure: Arb/FLINT (ball arithmetic), PARI/GP (`lfun`), eclib/mwrank (ranks, regulators, saturation), mpmath (probe work).
- Built in the Forge agentic workspace with continuous human direction, 2026-11-10.

## License

MIT — see [LICENSE](LICENSE).
