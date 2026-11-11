# M2 — Detector battery: the evaluator becomes an auditor

*2026-11-11 · v0.2.0 · campaign milestone M2*

M1 made a single certified value you could trust. M2 asks the harder operational question: **run the whole loop per curve, stamp every number with provenance, and let independent detectors demand agreement before anything is believed.**

## The pipeline (`scripts/detector.py`)

For one curve, end to end, every stage recorded with exact command, exit code, elapsed time, and sha256 of every input/output file:

1. **aₙ stream** — `gp`/`ellan`, Hecke-checked by the evaluator's gates later; stream length scales with conductor (v7 formula) so the selftests converge at any N ≤ 500,000.
2. **mwrank (eclib)** — unconditional rank proof + regulator on the saturated basis.
3. **PARI** — root number, Tamagawa product (via `elllocalred` ×2 when Δ>0), torsion, `ellrank` bounds.
4. **lane_b_cert v7** — certified [L(1), L′(1), L″(1)] balls + law/Hecke selftests (exit 4/3 on violation).
5. **gp `lfun` refs at `\p digits+300`** — fresh, per-run.
6. **Containment gate** — every reference must lie inside its certified ball (M1 rule: refs out-precision balls by hundreds of digits).

Then the detectors, ordered by cost:

- **D1** — certified non-vanishing: the L^(k) ball (k matched to root number) excludes 0 ⇒ ord ≤ k.
- **D3** — parity: provable algebraic rank parity must match the root-number parity.
- **D2** — implied-|Sha| square test (rank 2, ε=+1): a certified interval for

  |Sha| = L″(1)·T² / (2!·Ω·Reg·∏c_p)

  must contain a perfect-square integer, else ALARM (contradicts BSD formula + Sha finite + analytic rank 2 — rule #2 then demands independent evidence before any alarm is believed).

## Validation batch (reproducible)

`python3 scripts/detector.py --sample 10 --seed 20261110 --digits 100` — seeded pick of 10 rank-2 curves from the LMFDB rank-2 export, plus the 389a1 calibration curve through the same pipeline. **11/11 CLEAN:**

- containment True on every curve
- mwrank rank 2 unconditional == PARI ellrank bounds [2,2] on every curve
- regulator cross-check (mwrank's saturated-basis regulator vs the height-matrix determinant at high precision): agreement to machine precision (≤ 2.5e−14 relative)
- D1 certifies ord ≤ 2 everywhere; D3 parity consistent everywhere
- **D2: implied |Sha| window contains exactly square [1] on all 11 curves**

Runtime: 2.0 s (N = 389) to ~55 s (N ≈ 500k) per full audit at d100.

## What the battery caught during this batch

**Two real defects, both fixed, both by detector-designed process:**

1. **lane_b_cert v6 hardcoded small conductors.** The law selftest summed a fixed M = 3,000 terms (converges only for N ≲ 1,700 — at u = 0.5 the per-term decay is e^(−πn/√N)), and `n_eval` used `min()` where the decay-aware branch was meant to be a floor (pinning every big-N curve to 944 terms); true series tails at N ≈ 10⁵–10⁶ would have been ~e−8, i.e. garbage-wide certified balls. Every gate fired **loudly and correctly refused to certify** (exit 4, exactly as designed — the perturbation-tested gates working on real input). **v7 fix:** adaptive `M = ceil((dps+46)·ln10·√N/π) + 64`; `n_eval = max(8·(digits+18), ceil((digits+18)·ln10·√N/(2π)))`. Both preserve every v0.1.0-configuration certificate **byte-identically** (389a1 d100: n_eval still 944, containment mid−ref reproduced exactly). Design note: the law gate must evaluate **both sides directly** — a modular-transform shortcut would make both sides trivially equal and the perturbation gate vacuous.

2. **D2 bookkeeping had T² in the denominator.** Rank-2 BSD reads L″(1)/2! = Ω·Reg·|Sha|·∏c_p/T², so |Sha|_implied = L″·T²/(2!·Ω·Reg·∏c_p). The first implementation divided by T² instead — invisible on every calibration curve (trivial torsion or rank ≠ 2), and it fired the battery's first two alarms on the first two torsion-2 curves of the batch (485982.c1, 153990.l1: exact-1/16 windows = true value 1 shifted by T⁴ = 16). Post-fix both read |Sha| = 1.000… (38 digits). The formula was independently anchored on 11a1's exact known numbers (L(11a1,1) = Ω·5/25, c₁₁ = 5, T = 5).

**The lesson is the battery working as intended:** every alarm got investigated with independent evidence (the 11a1 anchor, the T⁴ rationality signature of the exact 1/16 windows) before any alarm was believed — and the alarm that looked like a finding was our own bookkeeping.

## Artifacts

- [certificates/m2/](../certificates/m2/) — 11 per-curve reports + `summary-seed20261110-d100.json` (regenerated post-fix; uniform provenance: detector sha `cfaafa13d67f1029`, lane_b_cert sha `1ba065c71c3cb07c`)
- The LMFDB rank-2 export (348,672 curves, conductors < 500,000; downloaded from lmfdb.org 2026-11-10) goes into `work/lmfdb-rank2-export.csv`; `mwrank` from an eclib build is expected at `~/.local/eclib/bin/mwrank`

## Next

- **Tier-0:** a 100-curve bit-agreement bridge (certified balls vs independent high-precision references, batched per-curve) — last gate before scale.
- **M3 (Tier-1):** the first verification tables over the full LMFDB rank-2 family: 348,672 curves / 279,056 isogeny classes / 129,892 conductors — measured ~2–55 s/audit at d100 ⇒ a few CPU-days sequential, parallelizable.
