# BSD Attack Plan — Certificate Shapes, Detectors, Milestones

Prize rules reminder (CMI, rev. 2018): publication in a Qualifying Outlet + **2 years** + general acceptance; either-direction results allowed for P vs NP and Navier–Stokes; for the others a counterexample is evaluated as a proposed solution; CMI decides in private. Solve today ⇒ earliest collection ~2028–2029.

## 1. Goal hierarchy (honest)

1. **PRIMARY:** a prize-grade counterexample certificate for BSD (rank part, Sha finiteness, or formula).
2. **SECONDARY (highly likely outcome):** certified verification tables for rank ≥ 2 curves — publishable mathematics regardless of the lottery.
3. **TERTIARY:** the certified toolchain itself (an interval-arithmetic BSD checker does not exist publicly in polished form — genuinely useful collateral).

## 2. Certificate shapes (what a win literally looks like)

Each is finite, machine-checkable, and refutes BSD-as-posed:

**CERT-1 (rank part, cheapest to check):** a curve `E/Q` with
- a displayed rational point `P` of infinite order (non-torsion certificate: reductions mod a finite set of primes with strictly increasing orders of `P mod p` — elementary and decisive), and
- a **certified** interval evaluation proving `L(E,1) ≠ 0` (⇒ `r_an = 0`).

Then `r_alg ≥ 1 ≠ 0 = r_an`: rank part of BSD false. Also contradicts the parity conjecture when `w(E) = -1` (expected for odd-rank curves), implying non-torsion Sha cascades. Prior probability per curve: astronomically low. Check cost per curve: trivial once the certified L-evaluator exists. Perfect swarm workload.

**CERT-2 (formula part, flagship):** a curve `E/Q` with
- `rank E(Q) = 2` **proven exactly** (2-descent upper bound + two displayed independent non-torsion points; regulators computed from certified generators),
- `Ω`, `tors`, `c_p` computed exactly (Tate's algorithm),
- `L''(E,1)/(2!)` evaluated with certified error bounds, and
- implied value `|Sha|_impl = |tors|² · L''(1)/(2!·Ω·Reg·∏c_p)` proven by intervals to be **not a perfect square integer**.

Since Cassels ⇒ finite Sha has square order, this refutes {BSD formula ∧ Sha finite}. This is the same certificate philosophy Boothby–Bradshaw ran in reverse (they verified to 10,000 digits *assuming Sha = 1*); we hunt curves where the assumption is impossible.

**CERT-3 (Sha finiteness failure, hardest):** proof that `Sha(E)[p^n]` is unbounded, or any rigorous demonstration of non-torsion Sha. No example is known for any curve over ℚ. Long-shot research lane, mostly theory; deprioritized but not zero (Selmer ratios over towers + Dokchitser-style Iwasawa tools).

**Anomalies (not certificates, still valuable):** extreme-precision near-misses — implied Sha that *should* be a square but is off by ~1e-40 (would indicate our own error first, then, if persistent across two implementations, a discovery-grade event to escalate).

## 3. The detector battery

| ID | Detector | Runs on | Certificate it can yield | Cost/curve |
| --- | --- | --- | --- | --- |
| D1 | rank ≥ 1 provable + certified `L(1) ≠ 0` | every rank ≥ 1 candidate | CERT-1 | very low |
| D2 | proven rank 2 + implied `\|Sha\|_impl` not a perfect square integer (certified) | rank-2 curves | CERT-2 | high (O(prec³)) |
| D3 | parity alarm: `w(E) = +1` with provable odd algebraic rank (or vice versa) | all curves | CERT-1-adjacent (refutes parity thms' hypotheses ⇒ Sha non-torsion) | very low |

D1/D3 are sieves (billions of cheap checks); D2 is the deep-dive lens. Every alarm must be re-derived by a second independent implementation (rule #2 in README).

## 4. Sieve tiers for candidate generation

- **Tier 0 — calibration:** known LMFDB rank-2 curves (e.g. 389a1) with published values. Purpose: prove our pipeline reproduces the literature exactly. Exit: bit-level agreement with Stein/LMFDB numbers on ≥ 100 curves.
- **Tier 1 — the complete map:** all rank ≥ 2 curves, conductors < 500,000 (complete LMFDB region — ~low-thousands of classes [EST], exact counts pulled in M0). Every one gets: exact rank proof, root number, certified L-values, D1+D3 first, D2 on rank-2 members. Expected outcome: a complete rigorous verification table for a region where currently *zero* curves have full verified BSD.
- **Tier 2 — prime-conductor wilderness:** conductor ≤ 300,000,000 prime-conductor curves (millions, mostly never rank-checked). Mestre–Nagao sums first (cheap pre-sieve for high-rank/high-Selmer suspects), then mwrank rank bounds, then D1/D3. High throughput, low per-curve depth.
- **Tier 3 — synthesis lane:** Mestre–Nagao-style parametric family sweeps for curves engineered toward high rank with even root number (maximizing CERT-1 prior against parity heuristics). Research-y; time-boxed.

## 5. Milestones

- **M0 — Toolchain + calibration.** Install PARI/GP, SageMath, eclib/mwrank, Arb bindings. Reproduce 389a1: rank 2 exactly, points `(0,0), (−1,1)`, `c_389 = 1`, regulator, `L''(1)/(2!·Ω)` to ≥ 1,000 digits, Boothby–Bradshaw 10,000-digit reproduction (budget ~1 CPU-week; start early, run in background). Also: pin down which three Millennium Problems OpenAI's 2026-10-06 release touched. Exit criteria: calibration table committed + docs updated with exact LMFDB rank-2 counts.
- **M1 — Certified L-derivative evaluator.** Interval/ball-arithmetic wrapper giving rigorous `L^(j)(E,1)` bounds. Validation: values agree with PARI/GP `lseries` (non-rigorous) at increasing precision, and error bounds hold under independent truncation. Exit: certified `L(1)`, `L'(1)`, `L''(1)` on the calibration set with proof-grade error bounds.
- **M2 — Detection battery, single-curve full loop.** D1, D3, D2 as runnable pipeline with audit trail (every number's provenance recorded). Exit: full report generated for 389a1 + 10 random rank-2 curves, independently reviewed.
- **M3 — Tier 1 campaign.** Complete scan of rank ≥ 2, conductors < 500,000. Exit: verification table published internally; zero unresolved alarms.
- **M4 — Tier 2 scale-out + alarm triage.** Swarm the prime-conductor wilderness. Exit: throughput SLO met and sustained.
- **M5 — Deep dives + disclosure gate.** Certified extreme-precision D2 on any survivor; two implementations agree before anything is called an alarm; no external disclosure or submission without explicit user sign-off (CMI clock only starts at publication).

## 6. Risks and walls (read before optimizing)

1. **Conjecture is probably true.** Delaunay heuristics: ~98%+ of rank ≥ 2 curves have trivial Sha; parity holds in every checked regime. Expected primary outcome: no counterexample. The verification tables are the real deliverable.
2. **The certification wall is real:** exact analytic rank beyond 0 (and, via Kolyvagin, 1) cannot be certified numerically in general — zeros can't be proven by interval containment. This is precisely why CERT-1/2 are designed around *non-vanishing* and *non-square* statements, the two shapes interval arithmetic CAN prove.
3. **Compute wall:** D2 precision costs O(prec³). 10,000 digits ≈ 1 CPU-week (Boothby–Bradshaw). Budget deep dives, don't spray them.
4. **Sha finiteness shadow:** even a "clean" implied-Sha integer of square shape for a rank-2 curve does not verify BSD (Sha finiteness unproven) — the tables must state conditional claims exactly. Precision in wording is part of the math.
5. **The OpenAI variable:** their 10-06 batch claimed progress on three unnamed Millennium Problems; if one touches BSD, M0 reconciles our plan with whatever is public before M3 begins.
6. **Priority hygiene (post-NS-dispute lesson):** timestamp every result, keep hashes in git, document what we derived vs. what we adopted from public sources — the Buckmaster/Alpöge dispute is the cautionary tale.

## 7. M0 kickoff — status (2026-11-10)

1. Toolchain: **PARI/GP 2.19.0 installed + verified** (Homebrew, arm64) — `lfun` derivative evaluation is the calibration workhorse. eclib: building from source (`work/eclib-src`; needs Homebrew prefix flags + CPATH for gmp headers — `make` loop closed, rebuild running). SageMath: brew cask **blocked** (needs interactive sudo) — alternatives recorded (micromamba userspace / source build), decision deferred to M2 prep.
2. LMFDB exact rank-2 counts: **blocked** (search-page count is JS-lazy + bot-gated; API challenged curl). Route recorded: pull the bulk export locally in M1 prep (needed for the M3 scan anyway).
3. 389a1 calibration: **PASSED** — see [results/m0-calibration.md](../results/m0-calibration.md). Key: Ω_BSD = 2·ω1_PARI (archimedean factor for two-component curves); ∏ c_p from elllocalred (NOT elltamagawa — PARI's global Tamagawa includes c_∞); |Sha|_impl = 1.000… exact to all 250 displayed digits.
4. OpenAI 3-problem identification: **DONE** — Riemann (003 quasi-RH), BSD (002 full formula, Selmer corank ≤ 1, incl. Sha finiteness), Hodge (032 CM abelian varieties + K3 products). P vs NP and Yang–Mills: untouched (overview section scans).
5. Precision pilot: **RUNNING** (`scripts/pilot2.sh` → `logs/precision-pilot.log`; 56 ms @ 300 digits, 33 ms @ 500 for L″(1); 10,000-digit Boothby–Bradshaw reproduction target — feasibility upgraded from CPU-weeks to likely hours/days).
