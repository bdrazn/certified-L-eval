#!/usr/bin/env python3
"""harness_agree.py — M1 cl-4: cross-lane agreement at 100 digits + perturbation gates.

Lane B = lane_b_cert.py certified balls (python-flint/Arb), invoked via the CLI.
Lane C = PARI/GP lfun values, recomputed fresh at digits+300 dps inside this run.
(d500 lesson: the ball radii (~e-584) can be TIGHTER than gp's own error at
\p digits+80 (~e-582) — refs must out-precision the ball radius by ~200 digits,
else containment fails spuriously. gp self-diff p580 vs p800 measured 5.6e-582.)

PASS per curve: the fresh lane-C value lies INSIDE the lane-B ball
   (|mid - ref| <= (hi-lo)/2 follows automatically from containment).
Perturbation: the sign-corrupted a_n stream must FAIL (nonzero exit:
Hecke gate or pointwise-law gate).

Usage: python3 scripts/harness_agree.py [--digits 100] [--curves all]
Record: logs/harness-d100.log, work/harness-d100.json
"""
import json, os, re, subprocess, sys, time

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
os.chdir(ROOT)

from mpmath import mp
mp.dps = 460
def mpv(b):
    return mp.mpf(str(b).strip("[] ").split("+/-")[0])

CURVES = [
    dict(tag="389a1", N=389, an="work/an-389a1.txt", eps=1,  k=2, prefix="1.5186", refkey="389_L2"),
    dict(tag="37a1",  N=37,  an="work/an-37a1.txt",  eps=-1, k=1, prefix="0.3059", refkey="37_L1"),
    dict(tag="20a1",  N=20,  an="work/an-20a1.txt",  eps=1,  k=0, prefix="0.4707", refkey="20_L0"),
]

GP = r"""
default(parisizemax, 4000000000);
default(threadsizemax, 1073741824);
\p 400
E = ellinit([0,1,1,-2,0]);
lf = lfuncreate(E);
print("389_L2 = ", Str(lfun(lf, 1, 2)));
E = ellinit([0,0,1,-1,0]);
lf = lfuncreate(E);
print("37_L1 = ", Str(lfun(lf, 1, 1)));
E = ellinit([0,1,0,-1,0]);
lf = lfuncreate(E);
print("20_L0 = ", Str(lfun(lf, 1)));
print("GPREF-DONE");
quit
"""

def gp_refs(prec=None):
    script = GP.replace("p 400", "p %d" % prec) if prec else GP
    p = subprocess.run(["gp", "-q"], input=script, capture_output=True, text=True, timeout=1800)
    refs = {}
    for line in p.stdout.splitlines():
        line = line.strip()
        m = re.match(r"(389_L2|37_L1|20_L0) = ([+\-0-9.Ee]+)$", line)
        if m:
            refs[m.group(1)] = m.group(2)
    assert "GPREF-DONE" in p.stdout and len(refs) == 3, "gp ref failed: %r" % p.stdout[-200:]
    return refs

def run_lane_b(tag, N, an, eps, k, digits, extra=()):
    out = "work/hb-%s-d%d.json" % (tag, digits)
    t0 = time.time()
    p = subprocess.run([sys.executable, "scripts/lane_b_cert.py", "--N", str(N),
                        "--an", an, "--eps", str(eps), "--k", str(k),
                        "--digits", str(digits), "--json-out", out, *extra],
                       capture_output=True, text=True, timeout=7200)
    return p, out, time.time() - t0

def main():
    digits = 100
    if "--digits" in sys.argv:
        digits = int(sys.argv[sys.argv.index("--digits") + 1])
    fresh = "--fresh" in sys.argv
    mp.dps = 3 * digits + 160
    ref_prec = digits + 300
    refs = gp_refs(prec=ref_prec)
    print("[gp refs at \\p %d]:" % ref_prec, {k: v[:44] for k, v in refs.items()})
    results, ok = [], True
    for c in CURVES:
        out = "work/hb-%s-d%d.json" % (c["tag"], digits)
        reused = os.path.exists(out) and not fresh
        if reused:
            print("[laneB-reuse] %s -> %s (certified ball on disk; --fresh recomputes)" % (c["tag"], out))
            dt = 0.0
        else:
            p, out, dt = run_lane_b(c["tag"], c["N"], c["an"], c["eps"], c["k"], digits)
            if p.returncode != 0:
                print("[FAIL-laneB] %s exit=%d %s" % (c["tag"], p.returncode, p.stderr[:300]))
                ok = False
                continue
        txt = open(out).read()
        j = json.loads(txt[txt.index("{"):txt.rindex("}") + 1])
        balls = re.findall(r"\[([+\-0-9.]+(?:e-[0-9]+)?) \+/- ([0-9.]+e-[0-9]+)\]", txt)
        cand = [b for b in balls if b[0].startswith(c["prefix"])]
        assert cand, "%s: no ball with prefix %r" % (c["tag"], c["prefix"])
        lo, hi = mp.mpf(cand[0][0]), mp.mpf(cand[-1][0])
        ref = mp.mpf(refs[c["refkey"]])
        inside = (lo <= ref <= hi)
        midoff = (lo + hi) / 2 - ref
        rad = (hi - lo) / 2
        passed = bool(inside)
        ok = ok and passed
        print("[%s-k%d d%d] %.1fs  inside=%s  radius=%s  mid-ref=%s"
              % (c["tag"], c["k"], digits, dt, inside, mp.nstr(rad, 6), mp.nstr(midoff, 6)))
        results.append(dict(tag=c["tag"], k=c["k"], digits=digits, seconds=round(dt, 2), reused=reused,
                            inside=inside, radius=str(mp.nstr(rad, 30)), mid_ref=str(mp.nstr(midoff, 30)),
                            ref=refs[c["refkey"]], lo=str(mp.nstr(lo, 120)), hi=str(mp.nstr(hi, 120))))
    # perturbation gates on all three curves
    for c in CURVES:
        p, _, dt = run_lane_b(c["tag"], c["N"], c["an"], c["eps"], c["k"], digits, extra=("--perturb",))
        caught = (p.returncode != 0)
        ok = ok and caught
        print("[perturb-%s] exit=%d caught_corruption=%s (%.1fs) %s"
              % (c["tag"], p.returncode, caught, dt, (p.stderr or "").strip()[:60]))
        results.append(dict(tag=c["tag"], perturb=True, exit=p.returncode, caught=caught))
    json.dump(dict(digits=digits, ok=ok, results=results),
              open("work/harness-d%d.json" % digits, "w"), indent=1)
    print("HARNESS-%s" % ("PASS" if ok else "FAIL"))
    return 0 if ok else 1

if __name__ == "__main__":
    sys.exit(main())
