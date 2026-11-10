# BSD: Computational State of the Art (as of 2026-11-10)

Every load-bearing claim is sourced. Confidence labels: [PROVEN] = theorem; [COMPUTED] = peer/verifier-checked computation; [EST] = estimate/heuristic; [NEWS] = reported event.

## 1. What BSD says

Let `E/Q` be an elliptic curve, `L(E,s)` its Hasse–Weil L-function.

**Rank part:** the Mordell–Weil rank `r_alg` equals the analytic rank `r_an = ord_{s=1} L(E,s)`.

**Refined formula (for r = r_an):**

```
L^(r)(E,1) / r!  =  Ω_E · Reg(E/Q) · |Sha(E/Q)| · ∏_p c_p  /  |E(Q)_tors|²
```

where `Ω_E` = real period of the minimal model, `Reg` = regulator (height-pairing determinant), `c_p` = Tamagawa numbers, and `Sha` = Tate–Shafarevich group. BSD asserts three things at once: (a) the rank identity, (b) finiteness of Sha, (c) the formula. Refuting **any one** of the three refutes BSD as posed.

**Square trick (our main detector lever):** if Sha is finite, the Cassels–Tate pairing forces `|Sha|` to be a **perfect square** [PROVEN]. So on any curve where everything except Sha is provably computed, the implied value `|Sha|_impl = |tors|² · L^(r)(1) / (r! · Ω · Reg · ∏ c_p)` must be a perfect square integer. A certified non-square (or non-integer) implied value refutes the conjunction {BSD formula, Sha finite} — which refutes BSD.

## 2. What is proven / known

- [PROVEN] **Analytic rank ≤ 1 ⇒ BSD essentially known.** Gross–Zagier + Kolyvagin (with modularity, Wiles et al.): `r_alg = r_an` when `r_an ≤ 1`, plus an algorithm (not necessarily practical) to verify the formula (Stein AMS talk slides, 2010).
- [COMPUTED] Stein (2009, Math. Comp.): BSD formula holds at primes ℓ for all non-CM curves with conductor ≤ 1000, rank ≤ 1, `E[ℓ]` irreducible, `ℓ ∤ c_p` — via 2/3-descent, Kolyvagin/Kato Euler systems, Sage.
- [COMPUTED] Stein–Wuthrich: exact Sha values for all rank ≤ 1 curves of conductor up to 130,000 (as of a 2012 note; "might be more by now" —Alex B., math.stackexchange #150065). Unconditional `Sha(E)[p]` computations for many higher-rank curves.
- [COMPUTED] The rank-2 curve **389a1** (`y(y+1) = x(x−1)(x+2)`): Boothby–Bradshaw verified the BSD formula to 10,000 decimal digits *assuming `Sha = 1`* (~1 week CPU; cost O(prec³) ⇒ a million digits ≈ a million weeks). Stein–Wuthrich: `Sha(389a1)[p] = 0` for all p < 2466 except possibly p = 107, 599, 1049 (non-ordinary residue classes).
- [PROVEN] Parity machinery: the parity conjecture `(-1)^rank = w(E/K)` holds conditionally (e.g. finiteness of 2- and 3-primary Sha of relevant twists — Dokchitser & Dokchitser, arXiv:0906.1815). Open in general over ℚ.
- [PROVEN] Mazur torsion theorem bounds `E(Q)_tors`; descent gives provable rank upper bounds (this is how exact ranks are certified).

## 3. What is NOT proven (the live frontier)

- **No rank ≥ 2 curve has the full BSD formula known** — for even 389a1, Sha finiteness itself is open (Stein 2010: "BSD formula not known for a single curve of rank ≥ 2").
- **No rank ≥ 2 curve over ℚ has Sha proven finite** (Alex B., math.SE 2012, consistent with the literature since).
- Exact algebraic rank is unconditionally proven only up to **22** (Samuel Li, Sep 2026 — Dujella rank-history page). Higher ranks are known only as lower bounds: current record **rank ≥ 32** (Elkies–Klagsbrun–Li, 2026; previous: 31 and 30 Alpöge–Howell 2026, 29 Elkies–Klagsbrun 2024, 28 Elkies 2006).
- Heuristics [EST]: Delaunay's Cohen–Lenstra-style heuristics predict ≈98%+ of rank ≥ 2 curves have trivial Sha — the needle must be thin if it exists at all.

**Consequence:** the hunt territory is precisely {curves with provable algebraic rank ≥ 2}. Rank ≤ 1 curves are covered by theorem; rank ≥ 2 is 100% open, unswarmed, and certificate-shaped.

## 4. Databases and scale

Per LMFDB completeness page (revised 2023-07): the ℚ-curve database holds **3,824,372 curves in 2,917,287 isogeny classes**, conductor ≤ 299,996,953, comprising:

- all curves with conductor < 500,000 (Cremona, eclib/modular symbols);
- all curves with 7-smooth conductor (Matschke–von Känel, Bennett–Gherga–Rechnitzer);
- all prime conductors ≤ 300,000,000 (Stein–Watkins, extended and completed by Bennett–Gherga–Rechnitzer).

Rank ≥ 2 curves up to conductor 500,000 number in the low thousands of isogeny classes [EST — extract exact counts from LMFDB in M0]. Prime-conductor territory adds millions more candidate curves, most never examined for rank.

## 5. Tools (all free unless noted)

- **PARI/GP** — `ellsearch`, `ellheegner`, Dokchitser `lseries` (high-precision L-derivative evaluation), `ellrank`.
- **SageMath** — descent (`TwoDescent`, `ThreeDescent`), modular symbols, `E.gens()`, `E.tamagawa_numbers()`, `E.period_lattice()`, ball arithmetic via Arb.
- **eclib / mwrank** (Cremona) — industrial-strength rank provers via 2-descent.
- **Magma** (licensed) — `RankBounds`, Heegner points; best-in-class but paid.
- **Arb / ball arithmetic** — the certification layer for interval-checked L-values.

## 6. 2026 context (why now)

- OpenAI claimed Navier–Stokes (2026-09-08, Science) and published 377 math results (2026-10-06/07, Fortune/WaPo), claiming "progress but not full solutions" on **three other Millennium Problems** (identity of the three not yet pinned down — verify in M0 whether BSD is among them and what they claim).
- ≥3 of those 377 were withdrawn after errors (Fortune). Formalization culture is now the credibility moat.
- Rank records jumped 28 → 29 → 30 → 31 → 32 between 2024 and 2026 — the numerical frontier is moving fast, partly with AI assistance. Note: **Levent Alpöge** (rank ≥ 30, 31 co-holder) is the same mathematician in the Navier–Stokes priority dispute with OpenAI.
- Community backlash (Association for Human Mathematics letter) means independent, carefully verified contributions have unusually high value right now — prize or no prize.

## 7. OpenAI 2026-10-06 catalog — entries bearing on Millennium problems

Primary source: OpenAI Research Catalog overview (github.com/openai/math, `overview.pdf`, October 6, 2026; 372 result families / 722 manuscripts). Entries cross-checked against the overview text pulled 2026-11-10:

- **002 — The full Birch–Swinnerton-Dyer formula in Selmer coranks zero and one.** Claims the full BSD leading-term formula for every E/Q whose full ℓ-power Selmer group has corank 0 or 1 for some prime ℓ, **including finiteness of Sha** for those curves. With entry 006 (Goldfeld's conjecture: analytic ranks 0 and 1 each have density 1/2 among quadratic twists of every E/Q), the catalog states this yields "full BSD for a density-one set of quadratic twists of every elliptic curve over ℚ." Status: claimed, unverified. [NEWS]
- **003 — The quasi-Riemann hypothesis.** Every Dirichlet L-function (including ζ) is zero-free in Re(s) > 7/8; alternate proof of zero-free Re(s) > 11/12; uniform exclusion of Landau–Siegel zeros. [NEWS; unpacked by human experts — Quanta 2026-10-19]
- **032 — The rational Hodge conjecture for CM abelian varieties and products of K3 surfaces.** Rational Hodge for every complex CM abelian variety (all dimensions/codimensions); via Milne's theorems ⟹ Tate conjecture for abelian varieties over finite fields + Hodge standard conjecture in every characteristic; algebraicity of Kuga–Satake correspondences for K3 products. Together with **001** (Milne's rationality conjecture for specialized Hodge classes), this is the catalog's Hodge-front entry. [NEWS]
- **004 — Hilbert's tenth problem over ℚ**: no algorithm decides whether an integer-coefficient polynomial in any number of variables has a rational zero. [NEWS]
- **030 — Modularity of elliptic curves over imaginary quadratic fields** (every such curve modular). [NEWS]

Section scans of the overview: Theoretical computer science = algorithms/lower bounds only — **no P vs NP entry** (closest: a permanent border-determinantal-complexity lower bound, entry ~108). Mathematical physics (entries 260–284 scanned) — **no Yang–Mills existence/mass-gap entry**. ⇒ The "progress on three other Millennium Prize problems" (Fortune, 2026-10-07) is exactly: **Riemann (003), BSD (002), Hodge (032/001)**.

**Implications for our campaign:** if 002 survives verification, the algebraic-Selmer-corank ≤ 1 region is *claimed closed* — including Sha finiteness there — pending community scrutiny (context: ≥3 catalog entries already withdrawn for errors within days of release). Our target territory {algebraic rank ≥ 2 / Selmer corank ≥ 2 curves} is strictly outside every claimed region and remains 100% open. Tiers and detectors unchanged; the M3 sweep is better motivated than before: no other claimed result touches our scan region.

## Sources

- Science, "How an AI math breakthrough ignited a controversy", 2026-09-08 — https://www.science.org/content/article/how-ai-math-breakthrough-ignited-controversy
- Clay Mathematics Institute, Navier–Stokes Announcement, 2026-09-11 — https://www.claymath.org/news/navier-stokes-announcement/
- CMI Rules for the Millennium Prize Problems (rev. 2018-09-26) — https://www.claymath.org/millennium-problems/rules/
- Quanta, "AI Proved the Quasi-Riemann Hypothesis. What Happens Now?", 2026-10-19 — https://www.quantamagazine.org/ai-proved-the-quasi-riemann-hypothesis-what-happens-now-20261019/
- Fortune, "OpenAI publishes solutions to more than 370 outstanding math challenges", 2026-10-07 — https://fortune.com/2026/10/07/openai-math-controversy-solutions-370-outstanding-challenges-published-criticisms-celebration/
- Washington Post, 2026-10-07 — https://www.washingtonpost.com/technology/2026/10/07/openai-releases-progress-more-than-300-math-research-problems/
- LMFDB completeness page (rev. 2023-07) — https://www.lmfdb.org/knowledge/show/rcs.cande.ec.q
- Dujella, "History of elliptic curves rank records" — https://dujella.github.io/tors/rankhist.html
- Stein, "The Birch and Swinnerton-Dyer Conjecture: A Template", AMS 2010-01-14 — https://www.wstein.org/talks/2010-01-14-ams/2010-01-14-ams_bsd/bsd.pdf
- Dokchitser & Dokchitser, "Root numbers and parity of ranks of elliptic curves" — https://arxiv.org/html/0906.1815v1
- math.stackexchange #150065 (Sha data availability; Alex B. answer) — https://math.stackexchange.com/questions/150065/
