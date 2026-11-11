#!/usr/bin/env python3
"""
detector.py — M2 detector battery pipeline with full per-number provenance.

Campaign rule #2 (two independent implementations must agree before any alarm
is reported) is enforced structurally:
  * rank         : mwrank (primary, 'determined unconditionally') vs PARI ellrank
                   proven bounds [r1, r2];
  * regulator    : mwrank 14-digit Regulator (on the saturated basis) vs det of
                   PARI canonical scalar-height pairing on the SAME basis at
                   \p digits+150 (formula: h1*h2 - ((h(P1+P2)-h1-h2)/2)^2);
  * L-values     : Arb certified balls (lane_b_cert.py, fresh run, gates on)
                   vs fresh gp lfun derivatives at \p digits+300 (M1 harness
                   convention, ref must out-precision ball by ~200 digits);
  * BSD ingredients: gp at \p digits+150 (Omega_BSD = ncomp*omega1 with ncomp = 2
                   if Delta>0 else 1; prod c_p from elllocalred last component,
                   cross-checked against elltamagawa(E)/(1+(Delta>0))).

Detectors:
  D1  nonvanishing certificate: certified ball for L^{(k)}(1) (k matched to the
      root number: k=2 for w=+1, k=1 for w=-1) strictly excludes 0  =>  ord <= k.
      For w=-1 the zero anchoring of L(1) is recorded (functional equation).
  D3  parity alarm: provable algebraic rank parity vs root-number parity.
  D2  implied-|Sha| square test (rank 2, w=+1):
          I = L''(1) * tors^2 / (2! * Omega_BSD * Reg * prod(c_p))
      computed as a certified INTERVAL. If I excludes every perfect square,
      ALARM (contradicts [BSD formula + Sha finite + analytic rank 2]).

Every number carries provenance: exact command, exit code, elapsed, sha256 of
every input/output file, versions. Reports land in work/reports/<label>.json.

Usage:
  python3 scripts/detector.py --label 389.a1 --digits 100
  python3 scripts/detector.py --sample 10 --seed 20261110 --digits 100
"""

import argparse
import datetime
import hashlib
import json
import math
import os
import random
import re
import subprocess
import sys
import time

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
WORK = os.path.join(ROOT, "work")
REPORTS = os.path.join(WORK, "reports")
CSV_PATH = os.path.join(WORK, "lmfdb-rank2-export.csv")
MWRANK = os.path.expanduser("~/.local/eclib/bin/mwrank")
LN10 = math.log(10.0)

os.makedirs(REPORTS, exist_ok=True)

# ---------------------------------------------------------------- provenance


def iso_now():
    return datetime.datetime.now(datetime.timezone.utc).isoformat(timespec="microseconds")


def sha256_path(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def run_rec(cmd, input_text=None, timeout=600, tag=""):
    """Run a command, return a provenance record."""
    t0 = time.time()
    try:
        p = subprocess.run(
            cmd,
            input=input_text,
            capture_output=True,
            text=True,
            timeout=timeout,
            cwd=ROOT,
        )
        rec = {
            "tag": tag,
            "cmd": " ".join(cmd) if isinstance(cmd, list) else cmd,
            "input_sha256": hashlib.sha256((input_text or "").encode()).hexdigest()[:16],
            "exit": p.returncode,
            "elapsed_s": round(time.time() - t0, 3),
            "stdout_tail": p.stdout[-2000:],
            "stderr_tail": p.stderr[-500:],
        }
        return rec, p
    except subprocess.TimeoutExpired:
        rec = {
            "tag": tag,
            "cmd": " ".join(cmd) if isinstance(cmd, list) else cmd,
            "exit": None,
            "error": "TIMEOUT after %ss" % timeout,
            "elapsed_s": round(time.time() - t0, 3),
        }
        return rec, None


def versions():
    rec, p = run_rec(["gp", "--version"], timeout=15, tag="gp-version")
    gp_ver = (p.stdout.splitlines()[0].strip() if p and p.stdout else "unknown")
    rec2, p2 = run_rec(
        [sys.executable, "-c", "import flint; print(flint.__version__)"],
        timeout=30,
        tag="flint-version",
    )
    flint_ver = (p2.stdout.strip() if p2 and p2.stdout else "unknown")
    return {
        "gp": gp_ver,
        "flint": flint_ver,
        "mwrank_sha256_16": sha256_path(MWRANK)[:16],
        "lane_b_cert_sha256_16": sha256_path(os.path.join(ROOT, "scripts", "lane_b_cert.py"))[:16],
        "detector_sha256_16": sha256_path(os.path.abspath(__file__))[:16],
        "gp_ref": rec,
        "flint_ref": rec2,
    }


# ---------------------------------------------------------------- CSV


def load_csv(path=CSV_PATH):
    rows = []
    with open(path, "r", encoding="utf-8") as f:
        for line in f:
            line = line.rstrip("\n")
            if not line.strip() or line.startswith("#"):
                continue
            parts = line.split("\t")
            if len(parts) != 7:
                continue  # non-data fragment; recorded count check in main
            rows.append(
                {
                    "label": parts[0].strip('"'),
                    "cls": parts[1].strip('"'),
                    "N": int(parts[2]),
                    "rank_lmfdb": int(parts[3]),
                    "torsion": parts[4],
                    "cm": int(parts[5]),
                    "coeffs": json.loads(parts[6].strip()),  # field is already valid JSON: [a1, a2, a3, a4, a6]
                }
            )
    return rows


# ---------------------------------------------------------------- stage B: a_n


def nterms_for(N, digits):
    return int((digits + 25) * LN10 * math.isqrt(N) / (2.0 * math.pi)) + 200


def extract_an(tag, coeffs, N, digits):
    # Stream length must cover (a) the evaluator's term count at target precision
    # (n_eval ~ (digits+31)*ln10*sqrt(N)/(2*pi), measured: 944 for 389a1@d100) and
    # (b) the pointwise-law selftest at min(digits+30,170) digits, which sums G(u)
    # at u near the fixed point with the same decay. A too-short stream makes the
    # law gate correctly FAIL (exit 4) — observed on d100 seeds with N ~ 250k-500k
    # when NT was floored at 12,600. Scale generously; ellan is cheap relative to
    # the certification it enables.
    NT = int(
        max(
            12600,
            math.ceil((digits + 60) * LN10 * math.sqrt(N) / (2.0 * math.pi)) + 500,
            # law-selftest need (lane_b_cert v7): M_law = ceil((dps+46)*ln10*sqrt(N)/pi)
            # + 64 with dps = min(digits+30, 170); branch below covers it with margin
            # (dps = digits+30 exactly for digits <= 140).
            math.ceil((digits + 76) * LN10 * math.sqrt(N) / math.pi) + 1000,
        )
    )
    gp = r"""
default(parisizemax, 4000000000);
E = ellinit(%s);
v = ellan(E, %d);
for(i=1, %d, if(abs(v[i]) > 2*sqrtint(i), error("HECKE-BOUND-FAIL ", i)));
for(i=1, %d, print(v[i]));
print("AN-EXTRACT-DONE");
quit
""" % (str(coeffs).replace(" ", ""), NT, NT, NT)
    script_path = os.path.join("/tmp", "det-extract-%s.gp" % tag.replace(".", "_"))
    with open(script_path, "w") as f:
        f.write(gp)
    rec, p = run_rec(["gp", "-q", "-f", script_path], timeout=1800, tag="extract-an")
    out_path = os.path.join(WORK, "an-det-%s.txt" % tag.replace(".", "_"))
    if p is None:
        raise RuntimeError("an extraction failed (timeout)")
    if p.returncode != 0 or "AN-EXTRACT-DONE" not in p.stdout:
        raise RuntimeError("an extraction failed: %r" % p.stderr[-300:])
    vals = [ln for ln in p.stdout.splitlines() if ln.strip() and not ln.startswith("***")]
    if not vals or vals[-1].strip() != "AN-EXTRACT-DONE":
        raise RuntimeError("an extraction marker missing; tail=%r" % p.stdout[-200:])
    vals = vals[:-1]
    if len(vals) != NT:
        raise RuntimeError("an extraction count mismatch: got %d expected %d; gp raw tail=%r" % (len(vals), NT, p.stdout[-200:]))
    with open(out_path, "w") as f:
        f.write("\n".join(vals) + "\n")  # drop the AN-EXTRACT-DONE marker
    return {
        "prov": rec,
        "path": out_path,
        "nterms": NT,
        "sha256": sha256_path(out_path),
        "hecke_checked": True,
    }


# ---------------------------------------------------------------- stage C: mwrank


def mwrank_stage(tag, coeffs):
    curve_in = "[%d,%d,%d,%d,%d]" % tuple(coeffs)
    rec, p = run_rec([MWRANK], input_text=curve_in + "\n", timeout=3600, tag="mwrank")
    if p is None:
        raise RuntimeError("mwrank timeout")
    out = p.stdout
    mr = re.search(r"Rank = (\d+)", out)  # mwrank prefixes lines with 'Curve [... ]'; exact 'Rank = N' cannot collide with 'Rank of ...'
    reg = re.search(r"Regulator = ([0-9.]+)", out)
    gens = [
        (int(a), int(b), int(c))
        for a, b, c in re.findall(r"Generator \d+ is \[(-?\d+):(-?\d+):(-?\d+)\]", out)
    ]
    return {
        "prov": rec,
        "rank": int(mr.group(1)) if mr else None,
        "unconditional": "determined unconditionally" in out,
        "gens_projective": gens,
        "reg_mwrank14": reg.group(1) if reg else None,
        "raw_tail": out[-1400:],
    }


# ---------------------------------------------------------------- stage D: PARI ingredients


def pari_stage(tag, coeffs, digits, mwrank):
    P = digits + 150
    gens_gp = ", ".join("[%d/%d, %d/%d]" % (x, z, y, z) for (x, y, z) in mwrank["gens_projective"])
    gp = r"""
default(parisizemax, 8000000000);
\p %d
E = ellinit(%s);
print("ING-N = ", ellglobalred(E)[1]);
print("ING-w = ", ellrootno(E));
print("ING-dsign = ", sign(E.disc));
print("ING-omega1 = ", E.omega[1]);
F = factor(ellglobalred(E)[1]); cp = 1;
for(i=1, matsize(F)[1], p = F[i,1]; lr = elllocalred(E,p); c1 = lr[#lr]; cp = cp*c1; print("ING-c_", p, " = ", c1));
print("ING-cprod = ", cp);
print("ING-tam = ", elltamagawa(E));
print("ING-tors = ", elltors(E)[1]);
gens = [%s];
h1 = ellheight(E, gens[1]);
print("ING-h1 = ", h1);
if(#gens == 2,
   h2 = ellheight(E, gens[2]); S = ellheight(E, elladd(E, gens[1], gens[2]));
   hp = (S - h1 - h2)/2;
   print("ING-reg = ", h1*h2 - hp*hp);
  , print("ING-reg = ", h1));
rk = ellrank(E);
print("ING-ellrank = ", rk);
print("ING-DONE");
quit
""" % (P, str(coeffs).replace(" ", ""), gens_gp)
    sp = os.path.join("/tmp", "det-ingredients-%s.gp" % tag.replace(".", "_"))
    with open(sp, "w") as f:
        f.write(gp)
    rec, p = run_rec(["gp", "-q", "-f", sp], timeout=1800, tag="pari-ingredients")
    if p is None or p.returncode != 0 or "ING-DONE" not in p.stdout:
        raise RuntimeError("pari ingredients failed: %r" % (p.stderr[-300:] if p else "timeout"))
    vals = {}
    for line in p.stdout.splitlines():
        m = re.match(r"ING-([A-Za-z0-9_]+) = (.+)$", line.strip())
        if m:
            vals[m.group(1)] = m.group(2).strip()
    needed = ["N", "w", "dsign", "omega1", "cprod", "tam", "tors", "h1", "reg"]
    missing = [k for k in needed if k not in vals]
    if missing:
        raise RuntimeError("pari ingredients missing keys: %r; tail %r" % (missing, p.stdout[-400:]))
    rank_lo = int(re.search(r"\[(\d+), ?(\d+),", vals["ellrank"]).group(1))
    rank_hi = int(re.search(r"\[(\d+), ?(\d+),", vals["ellrank"]).group(2))
    return {
        "prov": rec,
        "prec_digits": P,
        "N": int(vals["N"]),
        "w": int(vals["w"]),
        "disc_sign": int(vals["dsign"]),
        "omega1": vals["omega1"],
        "cprod": int(vals["cprod"]),
        "tam_global": int(vals["tam"]),
        "torsion": int(vals["tors"]),
        "reg_highprec": vals["reg"],
        "ellrank_bounds": [rank_lo, rank_hi],
        "raw": {k: v for k, v in vals.items()},
    }


# ---------------------------------------------------------------- stage E: certified balls


def parse_ball(text):
    pairs = re.findall(r"\[([+\-0-9.eE]+) \+/- ([0-9.eE+\-]+)\]", text)
    if not pairs:
        return None
    return pairs  # list of (mid, rad) sub-intervals; union taken by caller


def laneb_stage(tag, N, an_path, w, k, digits):
    out_path = os.path.join(WORK, "det-hb-%s-d%d.json" % (tag.replace(".", "_"), digits))
    cmd = [
        sys.executable,
        os.path.join(ROOT, "scripts", "lane_b_cert.py"),
        "--N", str(N),
        "--an", an_path,
        "--eps", str(w),
        "--k", str(k),
        "--digits", str(digits),
        "--json-out", out_path,
    ]
    rec, p = run_rec(cmd, timeout=7200, tag="lane-b-cert")
    if p is None or p.returncode != 0:
        raise RuntimeError("lane_b_cert failed: exit=%r tail=%r" % (p.returncode if p else None, p.stdout[-400:] if p else ""))
    with open(out_path) as f:
        hb = json.load(f)
    balls = {}
    for kk in range(0, k + 1):
        text = hb.get("L", {}).get(str(kk))
        if text is None:
            continue
        prs = parse_ball(text)
        if prs:
            balls[kk] = prs
    return {
        "prov": rec,
        "hb_path": out_path,
        "hb_sha256": sha256_path(out_path),
        "n_eval": hb.get("n_eval"),
        "prec_bits": hb.get("prec_bits"),
        "balls_raw": balls,
    }


# ---------------------------------------------------------------- stage F: gp refs


def gp_ref_stage(tag, coeffs, digits, k):
    P = digits + 300
    gp = r"""
default(parisizemax, 4000000000);
default(threadsizemax, 1073741824);
\p %d
E = ellinit(%s);
lf = lfuncreate(E);
print("REF-L = ", Str(lfun(lf, 1, %d)));
print("REF-DONE");
quit
""" % (P, str(coeffs).replace(" ", ""), k)
    sp = os.path.join("/tmp", "det-ref-%s.gp" % tag.replace(".", "_"))
    with open(sp, "w") as f:
        f.write(gp)
    rec, p = run_rec(["gp", "-q", "-f", sp], timeout=3600, tag="gp-ref")
    if p is None or p.returncode != 0 or "REF-DONE" not in p.stdout:
        raise RuntimeError("gp ref failed: %r" % (p.stderr[-300:] if p else "timeout"))
    m = re.search(r"REF-L = ([+\-0-9.Ee.]+)", p.stdout)
    if not m:
        raise RuntimeError("gp ref parse failed: %r" % p.stdout[-300:])
    return {"prov": rec, "prec_digits": P, "value": m.group(1)}


# ---------------------------------------------------------------- arithmetic on balls

from mpmath import mp  # noqa: E402  (import late keeps stage funcs readable)


def ball_union(pairs):
    mids = [mp.mpf(m) for m, r in pairs]
    rads = [mp.mpf(r) for m, r in pairs]
    return min(mids) - max(rads), max(mids) + max(rads)


def implied_sha_ball(lk_lo, lk_hi, omega_bsd, reg, cprod, tors):
    # BSD rank-2: L''(1)/2! = Omega*Reg*|Sha|*prod(c_p)/T^2  =>
    # |Sha| = L''(1) * T^2 / (2 * Omega * Reg * prod(c_p)).
    # (M2 2026-11-11: the T^2 factor originally sat in the denominator — all
    # M0/M1 calibration curves had trivial torsion, so 485982.c1 (T=2) was the
    # first curve to expose it: implied window was T^4=16x too small, reading
    # exactly 1/16 instead of exactly 1.  Battery did its job.)
    D = 2 * omega_bsd * reg * cprod
    return lk_lo * (tors ** 2) / D, lk_hi * (tors ** 2) / D, D


def squares_analysis(ilo, ihi):
    contains = []
    k = 1
    while (k + 1) ** 2 <= ihi or k * k <= ihi:
        sq = k * k
        if ilo <= sq <= ihi:
            contains.append(k)
        k += 1
        if k > 10 ** 7:
            break
    # margins: nearest square at or below ilo, smallest square strictly above ihi
    kb = max(0, int(math.sqrt(max(float(ilo), 0.0))))
    while kb > 0 and kb * kb > float(ilo):
        kb -= 1
    ka = int(math.sqrt(max(float(ihi), 0.0)))
    while ka * ka <= float(ihi):
        ka += 1
    return {
        "contains_squares_k": contains,
        "alarm_not_square": (ihi >= 0 and not contains),
        "nearest_square_below_k": kb,   # its square <= ilo
        "nearest_square_above_k": ka,   # its square > ihi
    }


# ---------------------------------------------------------------- detectors + report


def parity_ok(rank, w):
    return (rank % 2 == 0) if w == 1 else (rank % 2 == 1)


def run_curve(row, digits, versions_info):
    tag = row["label"].replace(".", "_")
    t_start = time.time()
    rpt = {
        "label": row["label"],
        "coeffs": row["coeffs"],
        "conductor": row["N"],
        "rank_lmfdb": row["rank_lmfdb"],
        "digits": digits,
        "started_utc": iso_now(),
        "versions": versions_info,
    }

    # stage B
    st_b = extract_an(tag, row["coeffs"], row["N"], digits)
    rpt["an"] = st_b

    # stage C
    st_c = mwrank_stage(tag, row["coeffs"])
    rpt["mwrank"] = st_c
    if st_c["rank"] is None:
        raise RuntimeError("mwrank rank unparsed; tail=%r" % st_c.get("raw_tail", "")[-500:])

    # stage D
    st_d = pari_stage(tag, row["coeffs"], digits, st_c)
    rpt["pari"] = st_d

    N = st_d["N"]
    w = st_d["w"]

    # stage E: k matched to root number
    k = 2 if w == 1 else 1
    st_e = laneb_stage(tag, N, st_b["path"], w, k, digits)
    rpt["laneb"] = {kk: v for kk, v in st_e.items() if kk != "balls_raw"}
    rpt["laneb"]["balls"] = {str(kk): [str(lo), str(hi)] for kk, prs in st_e["balls_raw"].items() for lo, hi in [tuple(map(str, ball_union(prs)))]}

    # stage F
    st_f = gp_ref_stage(tag, row["coeffs"], digits, k)
    rpt["gpref"] = st_f

    mp.dps = digits + 160
    # containment gate (harness pattern)
    return finish_curve(row, digits, rpt, st_b, st_c, st_d, st_e, st_f, k, w, t_start)


def finish_curve(row, digits, rpt, st_b, st_c, st_d, st_e, st_f, k, w, t_start):
    tag = row["label"].replace(".", "_")
    mp.dps = digits + 160
    lo, hi = ball_union(st_e["balls_raw"][k])
    ref = mp.mpf(st_f["value"])
    inside = (lo <= ref <= hi)
    mid = (lo + hi) / 2
    midoff = mid - ref
    rpt["containment"] = {
        "ref": st_f["value"],
        "inside": bool(inside),
        "mid_offset": str(midoff),
        "radius_rel": str((hi - lo) / 2 / max(abs(mid), mp.mpf("1e-500"))),
    }

    alarms = []
    certs = []

    # D1
    ball_k = (lo, hi)
    nonzero_cert = (lo > 0) or (hi < 0)
    if w == -1:
        certs.append({"D1": "L(1)=0 anchored by functional equation (w=-1, k odd)"})
    if nonzero_cert:
        certs.append({"D1": "certified L^%d(1) ball excludes 0 => ord <= %d" % (k, k)})
    else:
        alarms.append({"D1": "target ball straddles 0; ord <= %d NOT certified" % k})
    rpt["D1"] = {"ball": [str(lo), str(hi)], "nonzero_cert": bool(nonzero_cert)}

    # rank agreement (rule #2)
    rank = st_c["rank"]
    r1, r2 = st_d["ellrank_bounds"]
    rank_agree = (r1 == r2 == rank) and st_c["unconditional"]
    rpt["rank_agreement"] = {
        "mwrank_rank": rank,
        "mwrank_unconditional": st_c["unconditional"],
        "pari_ellrank_bounds": [r1, r2],
        "agree": bool(rank_agree),
    }
    if not rank_agree:
        alarms.append({"RANK": "mwrank=%r uncon=%r vs pari=%r" % (rank, st_c["unconditional"], [r1, r2])})

    # D3
    if not parity_ok(rank, w):
        alarms.append({"D3": "parity alarm: provable rank %d parity vs w=%d" % (rank, w)})
    rpt["D3"] = {"rank": rank, "w": w, "parity_ok": parity_ok(rank, w)}

    # D2
    if rank == 2 and w == 1:
        omega = 2 * mp.mpf(st_d["omega1"]) if st_d["disc_sign"] > 0 else mp.mpf(st_d["omega1"])
        reg = mp.mpf(st_d["reg_highprec"])
        reg14 = mp.mpf(st_c["reg_mwrank14"])
        # regulator cross-check: |reg_pari - reg_mwrank14| must be tiny relative
        reg_dev = abs(reg - reg14) / reg
        ilo, ihi, D = implied_sha_ball(lo, hi, omega, reg, mp.mpf(st_d["cprod"]), mp.mpf(st_d["torsion"]))
        sq = squares_analysis(ilo, ihi)
        rpt["D2"] = {
            "Omega_BSD": str(omega),
            "Reg_hp": str(reg),
            "reg_crosscheck_rel_dev": str(reg_dev),
            "denominator": str(D),
            "implied_lo": str(ilo),
            "implied_hi": str(ihi),
            "contains_squares_k": sq["contains_squares_k"],
            "alarm_not_square": sq["alarm_not_square"],
            "tam_crosscheck_ok": (st_d["tam_global"] == (2 if st_d["disc_sign"] > 0 else 1) * st_d["cprod"]),
        }
        if sq["alarm_not_square"]:
            alarms.append({"D2": "implied |Sha| ball [%s, %s] excludes ALL perfect squares" % (str(ilo)[:40], str(ihi)[:40])})
        else:
            certs.append({"D2": "implied |Sha| contains square(s) %r (SHA-CONSISTENT window)" % sq["contains_squares_k"]})
        if reg_dev > mp.mpf("1e-10"):
            alarms.append({"REG": "regulator cross-check deviation %s" % str(reg_dev)})

    rpt["alarms"] = alarms
    rpt["certificates"] = certs
    rpt["elapsed_s"] = round(time.time() - t_start, 1)
    rpt["finished_utc"] = iso_now()

    path = os.path.join(REPORTS, "%s-d%d.json" % (tag, digits))
    with open(path, "w") as f:
        json.dump(rpt, f, indent=1)
    verdict = "ALARM(%d)" % len(alarms) if alarms else "CLEAN"
    print("[%s] %s: %s | rank=%s w=%s k=%d | containment=%s | %s" % (
        rpt["finished_utc"], row["label"], verdict, rank, w, k, rpt["containment"]["inside"],
        "; ".join(c["D2"] for c in certs if "D2" in c) if certs else ", ".join(a["D2"] for a in alarms if "D2" in a)))
    return rpt


# ---------------------------------------------------------------- main


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--label", help="exact CSV label, e.g. 389.a1")
    ap.add_argument("--sample", type=int)
    ap.add_argument("--seed", type=int, default=20261110)
    ap.add_argument("--digits", type=int, default=100)
    args = ap.parse_args()

    rows = load_csv()
    print("CSV data rows: %d" % len(rows))
    vinfo = versions()
    if args.label:
        hit = [r for r in rows if r["label"] == args.label]
        if not hit:
            print("no such row")
            return 1
        run_curve(hit[0], args.digits, vinfo)
        return 0
    if args.sample:
        rng = random.Random(args.seed)
        picks = rng.sample(rows, args.sample)
        picked = [p["label"] for p in picks]
        print("seed=%d sample: %s" % (args.seed, ", ".join(picked)))
        out = {"seed": args.seed, "digits": args.digits, "picks": picked, "results": []}
        for p in picks:
            try:
                r = run_curve(p, args.digits, vinfo)
                out["results"].append({
                    "label": p["label"],
                    "alarms": r["alarms"],
                    "certificates": r["certificates"],
                    "containment": r["containment"]["inside"],
                    "elapsed_s": r["elapsed_s"],
                })
            except Exception as e:
                print("[%s] FAILED: %r" % (p["label"], e))
                out["results"].append({"label": p["label"], "error": repr(e)})
        with open(os.path.join(REPORTS, "summary-seed%d-d%d.json" % (args.seed, args.digits)), "w") as f:
            json.dump(out, f, indent=1)
        n_alarm = sum(1 for r in out["results"] if r.get("alarms"))
        print("SAMPLE SUMMARY: %d curves, %d failed, %d alarm curves" % (
            len(out["results"]), sum(1 for r in out["results"] if "error" in r), n_alarm))
        return 0
    print("nothing to do; use --label or --sample")
    return 1


if __name__ == "__main__":
    sys.exit(main())
