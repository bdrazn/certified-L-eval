#!/usr/bin/env python3
"""Tier-1 full-family scan driver (M3) -- parallel, resumable, compact rows.

Runs the full per-curve audit pipeline (scripts/detector.py at digits=100)
over every curve of the LMFDB rank-2 export, exactly like the Tier-0 bridge
but at family scale.

Design:
- Queue built once: dedup labels from detector.load_csv(), precompute the
  exact CM order discriminant (class-number-1 j-invariant bijection) and the
  torsion order so workers never need coeffs in memory, sort by (conductor,
  label) ASCENDING so any stopping point leaves a completed conductor tier.
- Static shards: worker i handles queue entries with idx % W == i (identical
  conductor distribution per worker; ascending order per shard; no locks).
- Per curve: detector.py as an isolated subprocess (cwd=ROOT, own process
  group, 900 s timeout, killpg, one retry) -- same invocation as the bridge.
- Compact rows only: JSONL rows carry every check verdict + provenance; full
  report files are deleted after extraction unless (row not ok -> alarms/ or
  failures/) or idx % 1000 == 0 (subsample/ spot-check archive).
- Supervisor: spawns workers detached (start_new_session, nice 10), polls
  every 60 s, relaunches dead-but-incomplete shards (bounded restarts), keeps
  a pid lock against double-start, writes shardNN.done sentinels and
  scan-complete.json at the end. Resume-from-manifest makes it reboot-safe.

Usage:
  python3 scripts/tier1_scan.py --build-queue
  python3 scripts/tier1_scan.py --supervise --workers 3 --limit 3 --smoke
  python3 scripts/tier1_scan.py --supervise --workers 8
  python3 scripts/tier1_scan.py --worker 0 --workers 8
  python3 scripts/tier1_scan.py --progress
  python3 scripts/tier1_scan.py --stop
"""
import argparse
import datetime
import json
import os
import shutil
import signal
import subprocess
import sys
import time

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(ROOT, "scripts"))
import detector  # noqa: E402
from bridge_tier0 import cm_disc_exact, torsion_order, ulabel  # noqa: E402

WORK1 = os.path.join(detector.WORK, "tier1")
QUEUE = os.path.join(WORK1, "queue.jsonl")
SUBSAMPLE_EVERY = 1000
TO_SEC = 900
TIERS = [25000, 50000, 100000, 200000, 300000, 400000, 500000]
FALLBACK_MEAN_S = 39.41  # Tier-0 bridge family mean


def rootdir(smoke):
    return os.path.join(detector.WORK, "tier1-smoke") if smoke else WORK1


def rows_path(root, shard):
    return os.path.join(root, "rows-shard%02d.jsonl" % shard)


def sentinel_path(root, shard):
    return os.path.join(root, "shard%02d.done" % shard)


def frontier_path(root, shard):
    return os.path.join(root, "frontier-shard%02d.json" % shard)


def worker_log(root, shard):
    return os.path.join(root, "logs", "worker%02d.log" % shard)


def safe(fn, default):
    try:
        return fn()
    except Exception:
        return default


# ---------------------------------------------------------------- queue ---

def build_queue():
    os.makedirs(WORK1, exist_ok=True)
    rows = detector.load_csv()
    seen, items = set(), []
    for r in rows:
        lab = r["label"]
        if lab in seen:
            continue
        seen.add(lab)
        items.append({
            "label": lab, "cls": r["cls"], "N": r["N"],
            "rank_lmfdb": r["rank_lmfdb"],
            "torsion_order": safe(lambda: torsion_order(r["torsion"]), -1),
            "cm": r["cm"],
            "cm_disc": safe(lambda: cm_disc_exact(r["coeffs"]), "ERROR"),
        })
    items.sort(key=lambda x: (x["N"], x["label"]))
    tmp = QUEUE + ".tmp"
    with open(tmp, "w", encoding="utf-8") as f:
        for i, it in enumerate(items):
            it["idx"] = i
            f.write(json.dumps(it, sort_keys=True) + "\n")
    os.replace(tmp, QUEUE)
    stats = ", ".join("N<%d: %d" % (t, sum(1 for it in items if it["N"] < t)) for t in TIERS)
    print("QUEUE BUILT: %d unique curves (from %d CSV rows) | %s" % (len(items), len(rows), stats),
          flush=True)


def load_queue():
    items = []
    with open(QUEUE, "r", encoding="utf-8") as f:
        for line in f:
            line = line.strip()
            if line:
                items.append(json.loads(line))
    return items


# ------------------------------------------------- per-curve audit row ----

def run_detector(label, digits, logdir, attempt, timeout=TO_SEC):
    """detector.py for one curve in its own process group; kill the whole
    group on timeout so a wedged mwrank grandchild cannot leak. Same shape as
    bridge_tier0.run_one (900 s here)."""
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


def extract_row(rpt, it, det_sha16):
    """Compact verdict row from a fresh report -- mirrors bridge assemble()."""
    with open(rpt, encoding="utf-8") as f:
        r = json.load(f)
    ra = r["rank_agreement"]
    bounds = ra.get("pari_ellrank_bounds") or []
    rank_ok = bool(ra.get("agree")) and ra.get("mwrank_rank") == 2 and bool(ra.get("mwrank_unconditional"))
    if len(bounds) == 2:
        rank_ok = rank_ok and bounds[0] <= 2 <= bounds[1]
    lmfdb_ok = r.get("rank_lmfdb") == 2 == it["rank_lmfdb"]
    t_pari = r["pari"].get("torsion")
    if isinstance(t_pari, str):
        t_pari = json.loads(t_pari)
    to = t_pari if isinstance(t_pari, int) else torsion_order(t_pari)
    torsion_ok = to == it["torsion_order"]
    d_cm = it["cm_disc"]
    cm_ok = ((it["cm"] == 0) == (d_cm in (None, "SINGULAR")) and
             (d_cm in (None, "SINGULAR") or d_cm == it["cm"]))
    cont_ok = r["containment"].get("inside") is True
    d2 = r["D2"]
    d2_ok = (d2.get("alarm_not_square") is False) and \
            (1 in [int(x) for x in (d2.get("contains_squares_k") or [])])
    d3_ok = r["D3"].get("parity_ok") is True
    noal = r.get("alarms") == []
    sha_ok = r["versions"].get("detector_sha256_16") == det_sha16
    checks = {"rank_ok": rank_ok, "uncond": bool(ra.get("mwrank_unconditional")),
              "lmfdb_ok": lmfdb_ok, "torsion_ok": torsion_ok, "cm_ok": cm_ok,
              "containment": cont_ok, "d2_square": d2_ok, "parity_ok": d3_ok,
              "no_alarms": noal, "sha_ok": sha_ok}
    row = dict(checks)
    row.update(ok=all(checks.values()),
               bad_checks=",".join(k for k, v in checks.items() if not v),
               reg_rel_dev=d2.get("reg_crosscheck_rel_dev"),
               contain_radius_rel=r["containment"].get("radius_rel"),
               detector_sha256_16=r["versions"].get("detector_sha256_16"),
               mwrank_rank=ra.get("mwrank_rank"))
    return row


def write_frontier(root, shard, done, last_N, ok_n, alarm_n, fail_n, ts):
    with open(frontier_path(root, shard), "w", encoding="utf-8") as f:
        json.dump({"shard": shard, "done": done, "frontier_N": last_N, "ok": ok_n,
                   "alarms": alarm_n, "fails": fail_n, "ts": ts}, f)


# ---------------------------------------------------------------- worker --

def run_worker(shard, W, digits, limit, smoke):
    root = rootdir(smoke)
    logdir = os.path.join(root, "logs")
    os.makedirs(logdir, exist_ok=True)
    os.makedirs(os.path.join(root, "alarms"), exist_ok=True)
    os.makedirs(os.path.join(root, "failures"), exist_ok=True)
    os.makedirs(os.path.join(root, "subsample"), exist_ok=True)
    items = load_queue()
    mine = [it for it in items if it["idx"] % W == shard]
    if limit:
        mine = mine[:limit]
    done = set()
    rp = rows_path(root, shard)
    if os.path.exists(rp):
        with open(rp, encoding="utf-8") as f:
            for line in f:
                if line.strip():
                    done.add(json.loads(line)["label"])
    det_sha16 = detector.sha256_path(os.path.join(ROOT, "scripts", "detector.py"))[:16]
    todo = [it for it in mine if it["label"] not in done]
    print("worker %02d start: mine=%d already-done=%d todo=%d root=%s sha=%s" % (
        shard, len(mine), len(done), len(todo), root, det_sha16), flush=True)
    t_start = time.time()
    ok_n = alarm_n = fail_n = 0
    cnt = 0
    for it in todo:
        label = it["label"]
        final = None
        for attempt in (1, 2):
            rc, wall, outp, errp = run_detector(label, digits, logdir, attempt)
            final = {"rc": rc, "wall_s": round(wall, 2), "attempt": attempt, "outp": outp, "errp": errp}
            if rc == 0:
                break
        row = {"label": label, "idx": it["idx"], "N": it["N"], "cls": it["cls"],
               "rank_lmfdb": it["rank_lmfdb"], "torsion_order": it["torsion_order"],
               "cm": it["cm"], "cm_disc": it["cm_disc"], "digits": digits,
               "rc": final["rc"], "attempt": final["attempt"], "wall_s": final["wall_s"],
               "ts": datetime.datetime.utcnow().isoformat() + "Z"}
        rpt = os.path.join(detector.REPORTS, "%s-d%d.json" % (ulabel(label), digits))
        if final["rc"] == 0:
            if os.path.exists(rpt):
                row.update(extract_row(rpt, it, det_sha16))
            else:
                row.update(ok=False, bad_checks="report_missing", rc=-3,
                           why="rc=0 but report missing: %s" % rpt)
        else:
            row.update(ok=False, bad_checks="detector_rc", rc=final["rc"],
                       why="detector failure rc=%s after retry (err: %s)" % (final["rc"], final["errp"]))
        alarm = (final["rc"] == 0) and row.get("ok") is False and row.get("bad_checks") != "report_missing"
        fail = final["rc"] != 0 or row.get("bad_checks") == "report_missing"
        keep = alarm or fail or (it["idx"] % SUBSAMPLE_EVERY == 0)
        if keep and os.path.exists(rpt):
            sub = "alarms" if alarm else ("failures" if fail else "subsample")
            shutil.move(rpt, os.path.join(root, sub, os.path.basename(rpt)))
        elif os.path.exists(rpt):
            os.remove(rpt)
        if not keep:
            for p_ in (final["outp"], final["errp"]):
                if p_ and os.path.exists(p_):
                    os.remove(p_)
        if row.get("ok") is True:
            for junk in (os.path.join(detector.WORK, "det-hb-%s-d%d.json" % (ulabel(label), digits)),
                         os.path.join(detector.WORK, "an-det-%s.txt" % ulabel(label))):
                if os.path.exists(junk):
                    try:
                        os.remove(junk)
                    except OSError:
                        pass
        with open(rp, "a", encoding="utf-8") as f:
            f.write(json.dumps(row, sort_keys=True) + "\n")
            f.flush()
            os.fsync(f.fileno())
        cnt += 1
        if row.get("ok"):
            ok_n += 1
        elif alarm:
            alarm_n += 1
            print("worker %02d ALARM %s N=%d bad=%s" % (shard, label, it["N"], row.get("bad_checks")), flush=True)
        else:
            fail_n += 1
            print("worker %02d FAIL %s rc=%s (kept logs %s)" % (shard, label, final["rc"], final["errp"]), flush=True)
        write_frontier(root, shard, len(done) + cnt, it["N"], ok_n, alarm_n, fail_n, row["ts"])
        print("worker %02d [%d/%d] %s N=%d wall=%.1fs ok=%s" % (
            shard, cnt, len(todo), label, it["N"], final["wall_s"], row.get("ok")), flush=True)
    with open(sentinel_path(root, shard), "w", encoding="utf-8") as f:
        json.dump({"shard": shard, "completed": len(done) + cnt, "mine": len(mine),
                   "ok": ok_n, "alarms": alarm_n, "fails": fail_n,
                   "finished_utc": datetime.datetime.utcnow().isoformat() + "Z",
                   "wall_s": round(time.time() - t_start, 1)}, f, indent=1)
    print("worker %02d DONE: %d curves in %.1fs (ok=%d alarms=%d fails=%d)" % (
        shard, len(done) + cnt, time.time() - t_start, ok_n, alarm_n, fail_n), flush=True)


# ------------------------------------------------------------ supervisor --

def spawn_worker(shard, W, digits, limit, smoke, root):
    cmd = ["nice", "-n", "10", sys.executable, os.path.abspath(__file__),
           "--worker", str(shard), "--workers", str(W), "--digits", str(digits)]
    if limit:
        cmd += ["--limit", str(limit)]
    if smoke:
        cmd += ["--smoke"]
    lf = open(worker_log(root, shard), "ab")
    return subprocess.Popen(cmd, cwd=ROOT, stdin=subprocess.DEVNULL,
                            stdout=lf, stderr=subprocess.STDOUT, start_new_session=True)


def supervise(args):
    W, digits, limit, smoke = args.workers, args.digits, args.limit, args.smoke
    root = rootdir(smoke)
    os.makedirs(os.path.join(root, "logs"), exist_ok=True)
    if not os.path.exists(QUEUE):
        print("SUPERVISOR: queue missing -- build first (--build-queue)", flush=True)
        return 2
    pf = os.path.join(root, "supervisor.pid")
    if os.path.exists(pf):
        try:
            old = int(open(pf).read().strip())
            os.kill(old, 0)
            print("SUPERVISOR: already running (pid %d) -- refusing double start" % old, flush=True)
            return 3
        except (ValueError, ProcessLookupError, PermissionError):
            pass
    with open(pf, "w", encoding="utf-8") as f:
        f.write(str(os.getpid()))
    items = load_queue()
    print("SUPERVISOR start pid=%d workers=%d digits=%d smoke=%s limit=%s family=%d" % (
        os.getpid(), W, digits, smoke, limit, len(items)), flush=True)
    procs, lives, restarts = {}, {}, {}
    for shard in range(W):
        if not os.path.exists(sentinel_path(root, shard)):
            procs[shard] = spawn_worker(shard, W, digits, limit, smoke, root)
            lives[shard] = 1
    t0 = time.time()
    while True:
        time.sleep(60)
        for shard in range(W):
            if os.path.exists(sentinel_path(root, shard)):
                continue
            p = procs.get(shard)
            if p is not None and p.poll() is None:
                continue
            if p is not None and p.poll() is not None and p.returncode != 0:
                print("SUPERVISOR %ds: shard %02d exited rc=%s" % (
                    round(time.time() - t0), shard, p.returncode), flush=True)
            n = restarts.get(shard, 0) + 1
            restarts[shard] = n
            if n > 40:
                if procs.get(shard):
                    procs.pop(shard)
                print("SUPERVISOR: shard %02d exceeded 40 restarts -- left down (check logs/worker%02d.log)" % (
                    shard, shard), flush=True)
                continue
            procs[shard] = spawn_worker(shard, W, digits, limit, smoke, root)
            lives[shard] = lives.get(shard, 0) + 1
            print("SUPERVISOR %ds: (re)spawned shard %02d (life #%d, restarts %d)" % (
                round(time.time() - t0), shard, lives[shard], n), flush=True)
        if all(os.path.exists(sentinel_path(root, s)) for s in range(W)) and \
                all(p.poll() is not None for p in procs.values()):
            print("SUPERVISOR: ALL SHARDS DONE -> scan complete", flush=True)
            with open(os.path.join(root, "scan-complete.json"), "w", encoding="utf-8") as f:
                json.dump({"completed_utc": datetime.datetime.utcnow().isoformat() + "Z",
                           "restarts": restarts, "lives": lives}, f, indent=1)
            try:
                os.remove(pf)
            except OSError:
                pass
            return 0


# -------------------------------------------------------------- progress --

def tail_rows(path, k=50):
    if not os.path.exists(path):
        return []
    with open(path, "rb") as f:
        f.seek(0, os.SEEK_END)
        size = f.tell()
        f.seek(max(0, size - 20000))
        data = f.read().decode("utf-8", "replace").splitlines()
    recs = []
    for line in data[-k:]:
        try:
            recs.append(json.loads(line))
        except Exception:
            pass
    return recs


def progress(args):
    smoke = args.smoke
    root = rootdir(smoke)
    if not os.path.exists(QUEUE):
        print("queue missing", flush=True)
        return 2
    items = load_queue()
    total = len(items)
    W = args.workers
    fstats = {}
    for s in range(64):
        fp = frontier_path(root, s)
        if not os.path.exists(fp):
            continue
        try:
            fstats[s] = json.load(open(fp))
        except Exception:
            pass
        if os.path.exists(sentinel_path(root, s)):
            fstats.setdefault(s, {"shard": s, "done": 0, "frontier_N": 0, "ok": 0, "alarms": 0, "fails": 0})
            fstats[s]["done_flag"] = True
    if not W:
        W = max([s + 1 for s in fstats] + [0]) or (8 if not smoke else 3)
    done_total = sum(v.get("done", 0) for v in fstats.values())
    done_total += sum(1 for v in fstats.values() if v.get("done_flag"))
    alarms = sum(v.get("alarms", 0) for v in fstats.values())
    fails = sum(v.get("fails", 0) for v in fstats.values())
    oks = sum(v.get("ok", 0) for v in fstats.values())
    walls = []
    for s in fstats:
        if fstats[s].get("done_flag"):
            continue
        walls += [r.get("wall_s", FALLBACK_MEAN_S) for r in tail_rows(rows_path(root, s))[-30:]]
    mean_recent = sum(walls) / len(walls) if walls else FALLBACK_MEAN_S
    remaining = max(0, total - done_total)
    eta_s = remaining * mean_recent / max(1, W)
    print("TIER1 PROGRESS: %d/%d rows (%.2f%%) | ok=%d alarms=%d fails=%d | mean-wall(recent)=%.1fs | W=%s | ETA ~ %.1f days" % (
        done_total, total, 100.0 * done_total / max(1, total), oks, alarms, fails,
        mean_recent, W, eta_s / 86400.0), flush=True)
    if smoke:
        items_by_idx = {it["idx"]: it for it in items}
    for s in sorted(fstats):
        v = fstats[s]
        print("  shard%02d: done=%s frontier_N=%s ok=%s alarms=%s fails=%s %s" % (
            s, v.get("done"), v.get("frontier_N"), v.get("ok"), v.get("alarms"), v.get("fails"),
            "[SENTINEL]" if v.get("done_flag") else ""), flush=True)
    for t in TIERS:
        if all(v.get("frontier_N", 0) >= t for v in fstats.values()) and fstats:
            print("  TIER COMPLETE: conductors < %d (all shards past)" % t, flush=True)
    return 0


# ------------------------------------------------------------------ stop --

def stop_scan(args):
    me = os.getpid()
    killed = 0
    for pat in ("tier1_scan.py --supervise", "tier1_scan.py --worker"):
        out = subprocess.run(["pgrep", "-f", pat], capture_output=True, text=True).stdout.split()
        for s in out:
            try:
                pid = int(s)
            except ValueError:
                continue
            if pid != me:
                try:
                    os.kill(pid, signal.SIGTERM)
                    killed += 1
                except ProcessLookupError:
                    pass
    print("STOP: SIGTERM sent to %d tier1 process(es)" % killed, flush=True)
    return 0


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--build-queue", action="store_true")
    ap.add_argument("--worker", type=int)
    ap.add_argument("--workers", type=int)
    ap.add_argument("--digits", type=int, default=100)
    ap.add_argument("--limit", type=int, help="worker: cap curves this invocation (smoke)")
    ap.add_argument("--supervise", action="store_true")
    ap.add_argument("--smoke", action="store_true", help="use work/tier1-smoke root")
    ap.add_argument("--progress", action="store_true")
    ap.add_argument("--stop", action="store_true")
    args = ap.parse_args()
    if args.build_queue:
        return build_queue() or 0
    if args.progress:
        return progress(args)
    if args.stop:
        return stop_scan(args)
    if args.supervise:
        return supervise(args)
    if args.worker is not None and args.workers:
        return run_worker(args.worker, args.workers, args.digits, args.limit, args.smoke) or 0
    print("nothing to do: pass --build-queue | --supervise --workers W | --worker i --workers W | --progress | --stop")
    return 2


if __name__ == "__main__":
    sys.exit(main())
