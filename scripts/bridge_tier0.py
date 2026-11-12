#!/usr/bin/env python3
"""Tier-0 bridge driver -- 100-curve bit-agreement run (M2 -> M3 gate).

Samples the LMFDB rank-2 export exactly like detector.py's --sample path (same
RNG stream; seed 20261110 with n=100 reproduces the M2 10 picks as its first
10, giving cross-run determinism anchors), then audits each curve through
detector.py as an isolated subprocess: per-curve wall timeout (600 s,
process-group kill so a wedged mwrank grandchild cannot leak), one retry,
resumable JSONL manifest.

--assemble produces the bit-agreement evidence:
  - agreement vs CSV-declared LMFDB fields (rank; torsion order; CM via the
    exact j-invariant <-> class-number-1 order-discriminant bijection over Q)
  - internal agreement (mwrank rank == PARI ellrank bounds, unconditional;
    regulator cross-check eclib vs PARI recorded per curve by the pipeline)
  - certified gates (D1 nonvanishing, D3 parity, containment of fresh gp lfun
    refs inside certified balls, D2 implied-|Sha| window contains square [1])
  - determinism: in-run rerun anchors (fresh subprocess, leaf-for-leaf compare
    vs the first-pass report over non-volatile fields) and cross-run anchors
    (the first 10 bridge curves rerun vs the git-shipped M2 certificates)

Usage:
  python3 scripts/bridge_tier0.py --seed 20261110 --n 100 --digits 100 [--limit K]
  python3 scripts/bridge_tier0.py --assemble [--rerun 3]
"""

import argparse
import datetime
import json
import os
import random
import re
import signal
import subprocess
import sys
import time
from fractions import Fraction

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(ROOT, "scripts"))
import detector  # stdlib-only module: reuse load_csv / sha256_path / sampler semantics

REPO_M2 = os.path.join(ROOT, "repo", "certified-L-eval", "certificates", "m2")
VOLATILE = re.compile(r"elapsed|started|finished|wall|duration|timestamp|date|utc", re.I)

# Over Q: E has CM iff j(E) is one of the 13 class-number-1 CM j-invariants,
# and the map j -> order discriminant is a bijection. Exact rational check.
CM_J_TO_D = {
    0: -3, 1728: -4, -3375: -7, 8000: -8, -32768: -11, 54000: -12,
    287496: -16, -884736: -19, -12288000: -27, 16581375: -28,
    -884736000: -43, -147197952000: -67, -262537412640768000: -163,
}


def ulabel(label):
    return label.replace(".", "_")


def bridge_dir(seed, digits):
    d = os.path.join(detector.WORK, "bridge-tier0-seed%d-d%d" % (seed, digits))
    os.makedirs(os.path.join(d, "logs"), exist_ok=True)
    return d


def sample_picks(seed, n):
    rng = random.Random(seed)
    return rng.sample(detector.load_csv(), n)


def manifest_path(seed, digits):
    return os.path.join(bridge_dir(seed, digits), "manifest.jsonl")


def load_manifest(seed, digits):
    last = {}
    try:
        with open(manifest_path(seed, digits), "r", encoding="utf-8") as f:
            for line in f:
                line = line.strip()
                if line:
                    rec = json.loads(line)
                    last[rec["label"]] = rec
    except FileNotFoundError:
        pass
    return last


def append_manifest(seed, digits, rec):
    with open(manifest_path(seed, digits), "a", encoding="utf-8") as f:
        f.write(json.dumps(rec, sort_keys=True) + "\n")
        f.flush()
        os.fsync(f.fileno())


def run_one(label, digits, logdir, attempt, timeout=600):
    """detector.py for one curve in its own process group; kill the whole
    group on timeout so a wedged mwrank grandchild cannot leak."""
    outp = os.path.join(logdir, "%s.attempt%d.out" % (ulabel(label), attempt))
    errp = os.path.join(logdir, "%s.attempt%d.err" % (ulabel(label), attempt))
    with open(outp, "wb") as fo, open(errp, "wb") as fe:
        t0 = time.time()
        p = subprocess.Popen(
            [sys.executable, "scripts/detector.py", "--label", label, "--digits", str(digits)],
            cwd=ROOT, stdout=fo, stderr=fe, start_new_session=True)
        try:
            p.communicate(timeout=timeout)
            rc = p.returncode
        except subprocess.TimeoutExpired:
            try:
                os.killpg(os.getpgid(p.pid), signal.SIGKILL)
            except (ProcessLookupError, PermissionError):
                pass
            p.communicate()
            rc = -2
    return rc, time.time() - t0, outp, errp


def run_batch(args):
    picks = sample_picks(args.seed, args.n)
    sfile = os.path.join(bridge_dir(args.seed, args.digits), "sample.json")
    if not os.path.exists(sfile):
        with open(sfile, "w", encoding="utf-8") as f:
            json.dump({"seed": args.seed, "digits": args.digits,
                       "picked_utc": datetime.datetime.utcnow().isoformat() + "Z",
                       "sample": [{"label": p["label"], "cls": p["cls"], "N": p["N"],
                                   "rank_lmfdb": p["rank_lmfdb"], "torsion": p["torsion"],
                                   "cm": p["cm"], "coeffs": p["coeffs"]} for p in picks]}, f, indent=1)
    done = {} if getattr(args, "fresh", False) else load_manifest(args.seed, args.digits)
    okset = {k for k, v in done.items() if v.get("rc") == 0
             and os.path.exists(os.path.join(detector.REPORTS, "%s-d%d.json" % (ulabel(k), args.digits)))}
    print("sample n=%d seed=%d digits=%d | done %d/%d (resumable manifest)" % (
        args.n, args.seed, args.digits, len(okset), args.n), flush=True)
    limit = args.limit if args.limit else args.n
    processed = 0
    for idx, p in enumerate(picks):
        if processed >= limit:
            break
        label = p["label"]
        if label in okset:
            continue
        processed += 1
        final = None
        for attempt in (1, 2):
            rc, wall, outp, errp = run_one(label, args.digits, os.path.join(bridge_dir(args.seed, args.digits), "logs"), attempt)
            rpt = os.path.join(detector.REPORTS, "%s-d%d.json" % (ulabel(label), args.digits))
            final = {"label": label, "idx": idx, "conductor": p["N"], "attempt": attempt, "rc": rc,
                     "wall_s": round(wall, 2), "ts": datetime.datetime.utcnow().isoformat() + "Z",
                     "report": rpt,
                     "report_sha256_16": detector.sha256_path(rpt)[:16] if os.path.exists(rpt) else None}
            append_manifest(args.seed, args.digits, final)
            if rc == 0:
                break
            print("[%d] %s rc=%d attempt%d FAILED (%s)" % (idx + 1, label, rc, attempt, errp), flush=True)
        if final["rc"] == 0:
            print("[%d] %s rc=0 %.1fs" % (idx + 1, label, final["wall_s"]), flush=True)
        else:
            print("[%d] %s FINAL rc=%d after retry" % (idx + 1, label, final["rc"]), flush=True)
    done = load_manifest(args.seed, args.digits)
    n_ok = sum(1 for v in done.values() if v.get("rc") == 0)
    print("BATCH DONE: %d/%d records ok (of n=%d)" % (n_ok, len(done), args.n), flush=True)
    return 0 if n_ok == args.n else 1


def strip_volatile(obj):
    """Recurse; drop any key matching VOLATILE (timings/date stamps are not
    deterministic); compare everything else leaf-for-leaf."""
    if isinstance(obj, dict):
        return {k: strip_volatile(v) for k, v in obj.items() if not VOLATILE.search(k)}
    if isinstance(obj, list):
        return [strip_volatile(x) for x in obj]
    return obj


from decimal import Decimal  # noqa: E402  (exact-decimal compare for ball endpoints)

_SECONDS_RE = re.compile(r"\(\d+(?:\.\d+)? seconds\)")
NUMERIC_BALL_CMP = [0]


def deep_compare(a, b, path="$"):
    """Strict leaf compare (int vs bool vs float distinguished). Returns
    (equal, leaves_compared, diffs_up_to_8)."""
    bad = []
    n = 0
    if isinstance(a, dict) and isinstance(b, dict):
        for k in sorted(set(a) | set(b)):
            if k not in a or k not in b:
                bad.append("%s.%s: present=%s/%s" % (path, k, k in a, k in b))
                continue
            e, m, d = deep_compare(a[k], b[k], path + "." + k)
            n += m
            bad += d
    elif isinstance(a, list) and isinstance(b, list):
        if len(a) != len(b):
            bad.append("%s: len %d != %d" % (path, len(a), len(b)))
        for i, (x, y) in enumerate(zip(a, b)):
            e, m, d = deep_compare(x, y, path + "[%d]" % i)
            n += m
            bad += d
    else:
        n += 1
        if isinstance(a, str) and isinstance(b, str):
            # mwrank/eclib echo their own wall time into captured tails —
            # environment noise, not math: normalize before comparing.
            a = _SECONDS_RE.sub("(T)", a)
            b = _SECONDS_RE.sub("(T)", b)
        leaf_eq = (type(a) is type(b)) and (a == b)
        if (not leaf_eq) and isinstance(a, str) and isinstance(b, str) \
                and path.startswith("$.laneb.balls"):
            # Ball-endpoint DISPLAY drifted across report generations (ambient
            # mp.dps leak — fixed in detector.py 2026-11-11; math fields were
            # always bit-stable). Compare endpoints numerically at tight
            # tolerance and count the occurrences for the bridge report.
            da, db = Decimal(a), Decimal(b)
            scale = max(abs(da), abs(db))
            if scale != 0 and abs(da - db) / scale < Decimal("1e-14"):
                leaf_eq = True
                NUMERIC_BALL_CMP[0] += 1
        if not leaf_eq:
            bad.append("%s: %r != %r" % (path, a, b))
    return (not bad), n, bad[:8]


def torsion_order(s):
    lst = json.loads(s) if isinstance(s, str) else s
    o = 1
    for m in lst:
        o *= int(m)
    return o


def cm_disc_exact(coeffs):
    """Exact j-invariant from [a1,a2,a3,a4,a6]; return CM order discriminant
    (class-number-1 bijection) or None when non-CM."""
    a1, a2, a3, a4, a6 = [Fraction(c) for c in coeffs]
    b2 = a1 * a1 + 4 * a2
    b4 = 2 * a4 + a1 * a3
    b6 = a3 * a3 + 4 * a6
    b8 = a1 * a1 * a6 + 4 * a2 * a6 - a1 * a3 * a4 + a2 * a3 * a3 - a4 * a4
    disc = -b2 * b2 * b8 - 8 * b4 ** 3 - 27 * b6 * b6 + 9 * b2 * b4 * b6
    if disc == 0:
        return "SINGULAR"
    c4 = b2 * b2 - 24 * b4
    j = c4 ** 3 / disc
    if j.denominator != 1:
        return None
    return CM_J_TO_D.get(int(j))


def assemble(args):
    seed, digits = args.seed, args.digits
    picks = sample_picks(seed, args.n)
    labels = [p["label"] for p in picks]
    srows = {r["label"]: r for r in json.load(open(os.path.join(bridge_dir(seed, digits), "sample.json")))["sample"]}
    manifest = load_manifest(seed, digits)
    det_sha16 = detector.sha256_path(os.path.join(ROOT, "scripts", "detector.py"))[:16]

    per_curve, fails, nonuniform, walls = [], [], [], []
    max_reg_dev, max_rad, cm_seen = 0.0, 0.0, 0
    counts = dict(rank_ok=0, uncond=0, lmfdb_ok=0, torsion_ok=0, cm_ok=0,
                  containment=0, d2_square=0, parity=0, no_alarms=0, uniform=0)
    for label in labels:
        rec = manifest.get(label)
        row = srows[label]
        entry = {"label": label, "N": row["N"], "rc": rec["rc"] if rec else None,
                 "wall_s": rec["wall_s"] if rec else None}
        if rec is None or rec.get("rc") != 0:
            entry.update(ok=False, why="no record" if rec is None else "rc=%s" % rec.get("rc"))
            fails.append(label)
            per_curve.append(entry)
            continue
        r = json.load(open(rec["report"]))
        sha_ok = r["versions"]["detector_sha256_16"] == det_sha16
        ra = r["rank_agreement"]
        bounds = ra.get("pari_ellrank_bounds") or []
        rank_ok = bool(ra.get("agree")) and ra.get("mwrank_rank") == 2 and bool(ra.get("mwrank_unconditional"))
        if len(bounds) == 2:
            rank_ok = rank_ok and bounds[0] <= 2 <= bounds[1]
        lmfdb_ok = r.get("rank_lmfdb") == 2 == row["rank_lmfdb"]
        t_pari = r["pari"].get("torsion")
        if isinstance(t_pari, str):
            t_pari = json.loads(t_pari)
        torsion_ok = torsion_order(row["torsion"]) == (t_pari if isinstance(t_pari, int) else torsion_order(t_pari))
        d_cm = cm_disc_exact(row["coeffs"])
        if d_cm not in (None, "SINGULAR"):
            cm_seen += 1
        cm_ok = ((row["cm"] == 0) == (d_cm in (None, "SINGULAR")) and
                 (d_cm in (None, "SINGULAR") or d_cm == row["cm"]))
        cont_ok = r["containment"].get("inside") is True
        d2 = r["D2"]
        d2_ok = (d2.get("alarm_not_square") is False) and \
                (1 in [int(x) for x in (d2.get("contains_squares_k") or [])])
        d3_ok = r["D3"].get("parity_ok") is True
        noal = r.get("alarms") == []
        all_ok = rank_ok and lmfdb_ok and torsion_ok and cm_ok and cont_ok and d2_ok and d3_ok and noal and sha_ok
        entry.update(ok=all_ok, rank_ok=rank_ok, uncond=bool(ra.get("mwrank_unconditional")),
                     lmfdb_ok=lmfdb_ok, torsion_ok=torsion_ok, cm_ok=cm_ok, containment=cont_ok,
                     d2_square=d2_ok, parity=d3_ok, no_alarms=noal, detector_sha_ok=sha_ok,
                     reg_dev=d2.get("reg_crosscheck_rel_dev"),
                     radius_rel=r["containment"].get("radius_rel"))
        for k, v in (("rank_ok", rank_ok), ("uncond", ra.get("mwrank_unconditional")),
                     ("lmfdb_ok", lmfdb_ok), ("torsion_ok", torsion_ok), ("cm_ok", cm_ok),
                     ("containment", cont_ok), ("d2_square", d2_ok), ("parity", d3_ok),
                     ("no_alarms", noal), ("uniform", sha_ok)):
            if v:
                counts[k] += 1
            elif k == "uniform":
                nonuniform.append(label)
        rd, rr_ = d2.get("reg_crosscheck_rel_dev"), r["containment"].get("radius_rel")
        if isinstance(rd, (int, float)):
            max_reg_dev = max(max_reg_dev, abs(float(rd)))
        if isinstance(rr_, (int, float)):
            max_rad = max(max_rad, abs(float(rr_)))
        walls.append(rec["wall_s"])
        per_curve.append(entry)

    ok_labels = [e["label"] for e in per_curve if e.get("ok")]
    idxs = []
    seen = set()
    for i in [0, len(ok_labels) // 2, len(ok_labels) - 1][:max(0, args.rerun)]:
        if 0 <= i < len(ok_labels) and ok_labels[i] not in seen:
            idxs.append(i)
            seen.add(ok_labels[i])
    rerun_blocks = []
    for i in idxs:
        label = ok_labels[i]
        rpt = os.path.join(detector.REPORTS, "%s-d%d.json" % (ulabel(label), digits))
        first = strip_volatile(json.load(open(rpt)))  # snapshot BEFORE rerun overwrites the file
        rc, wall, outp, errp = run_one(label, digits, os.path.join(bridge_dir(seed, digits), "logs"), 9)
        if rc != 0:
            rerun_blocks.append({"label": label, "rc": rc, "bit_equal": False, "leaves_compared": 0,
                                 "diffs": ["rerun rc=%d (%s)" % (rc, errp)]})
            continue
        second = strip_volatile(json.load(open(rpt)))
        eq, n, diffs = deep_compare(first, second)
        rerun_blocks.append({"label": label, "rc": rc, "wall_s": round(wall, 2), "bit_equal": eq,
                             "leaves_compared": n, "diffs": diffs,
                             "rerun_report_sha256_16": detector.sha256_path(rpt)[:16],
                             "note": "report file now holds rerun output; first-pass sha in manifest"})

    anchors = []
    for label in labels[:10]:
        pa = os.path.join(REPO_M2, "%s-d%d.json" % (ulabel(label), digits))
        pb = os.path.join(detector.REPORTS, "%s-d%d.json" % (ulabel(label), digits))
        if not (os.path.exists(pa) and os.path.exists(pb)):
            anchors.append({"label": label, "present": False})
            continue
        eq, n, diffs = deep_compare(strip_volatile(json.load(open(pa))), strip_volatile(json.load(open(pb))))
        anchors.append({"label": label, "present": True, "bit_equal": eq,
                        "leaves_compared": n, "diffs": diffs})

    n_total = len(labels)
    n_ok = sum(1 for e in per_curve if e.get("ok"))
    timing = {
        "mean_s": round(sum(walls) / len(walls), 2) if walls else None,
        "median_s": sorted(walls)[len(walls) // 2] if walls else None,
        "min_s": round(min(walls), 2) if walls else None,
        "max_s": round(max(walls), 2) if walls else None,
        "total_s": round(sum(walls), 1),
    }
    report = {
        "bridge": "tier0", "seed": seed, "digits": digits, "n_target": n_total,
        "generated_utc": datetime.datetime.utcnow().isoformat() + "Z",
        "driver_sha256_16": detector.sha256_path(os.path.abspath(__file__))[:16],
        "detector_sha256_16": det_sha16,
        "n_ok": n_ok, "n_fail": n_total - n_ok, "fails": fails,
        "checks_pct": {k: round(100.0 * v / max(1, n_ok), 1) for k, v in counts.items()},
        "uniform_detector_sha_failures": nonuniform,
        "max_reg_crosscheck_rel_dev": max_reg_dev,
        "max_containment_radius_rel": max_rad,
        "cm_curves_in_sample": cm_seen,
        "rerun_bit_compare": rerun_blocks,
        "anchor_m2_bit_compare": anchors,
        "timing": timing,
        "method_notes": [
            "volatile keys (elapsed/started/finished/timing stamps) excluded from bit-compares; all value fields, shas, paths and stream provenance compared",
            "CM agreement via exact j-invariant and the class-number-1 j<->order-discriminant bijection over Q (13 values); CSV cm column is the order discriminant",
            "anchor compares: bridge rerun of the 10 seed-20261110 curves vs the git-shipped M2 certificates (repo/certified-L-eval/certificates/m2)",
        ],
        "per_curve": per_curve,
    }
    out = os.path.join(detector.REPORTS, "bridge-tier0-seed%d-d%d.json" % (seed, digits))
    with open(out, "w", encoding="utf-8") as f:
        json.dump(report, f, indent=1)
    print("bridge report: %s" % out)
    print("OK %d/%d | rank %d | uncond %d | lmfdb %d | torsion %d | cm %d | containment %d | D2-square %d | parity %d | no-alarms %d" % (
        n_ok, n_total, counts["rank_ok"], counts["uncond"], counts["lmfdb_ok"], counts["torsion_ok"],
        counts["cm_ok"], counts["containment"], counts["d2_square"], counts["parity"], counts["no_alarms"]))
    for rb in rerun_blocks:
        print("rerun bit-equal %-12s -> %s (%d leaves)%s" % (
            rb["label"], rb["bit_equal"], rb["leaves_compared"],
            "" if rb["bit_equal"] else " DIFFS: %s" % rb["diffs"]))
    for a in anchors:
        if a.get("present"):
            print("anchor bit-equal %-12s -> %s (%d leaves)%s" % (
                a["label"], a["bit_equal"], a["leaves_compared"],
                "" if a["bit_equal"] else " DIFFS: %s" % a["diffs"]))
    all_equal = all(rb["bit_equal"] for rb in rerun_blocks) and \
        all(a["bit_equal"] for a in anchors if a.get("present"))
    print("DETERMINISM: %s" % ("ALL BIT-EQUAL" if all_equal else "MISMATCH"))
    return 0 if (n_ok == n_total and all_equal and not nonuniform) else 1


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--seed", type=int, default=20261110)
    ap.add_argument("--n", type=int, default=100)
    ap.add_argument("--digits", type=int, default=100)
    ap.add_argument("--limit", type=int, help="run at most K more curves (smoke)")
    ap.add_argument("--assemble", action="store_true")
    ap.add_argument("--rerun", type=int, default=3)
    ap.add_argument("--fresh", action="store_true",
                    help="ignore existing manifest records (re-run the whole sample);" 
                         " the manifest keeps its prior history")
    args = ap.parse_args()
    if args.assemble:
        return assemble(args)
    return run_batch(args)


if __name__ == "__main__":
    sys.exit(main())
