#!/usr/bin/env python
"""A/B timing of two variants in one session, for a change measured against the commit it builds on.

    bench/run ab --a tip --b serial [--datasets H1_...,H3_...] [--threads 16,64] [--reps 3]

Runs alternate (a then b on even repetitions, b then a on odd ones) after the load settles, so drift within the
session affects both sides alike. Every pair must produce the same output hashes. Records go to
<work>/ab/<a>_vs_<b>.jsonl (rerunning skips finished runs); the summary gives median wall times, their ratio
and the median of each `[phase]` marker interval.
"""
import argparse
import json
import os
import statistics
import sys
import time
import types

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import bench  # noqa: E402

DEFAULT_DATASETS = "H1_bam_reps_karpas422,H3_bam_ct_gm12878_grocap,T4_bam_reps_hct116,T5_bam_ct_k562_grocap"


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--a", required=True, help="baseline variant")
    ap.add_argument("--b", required=True, help="variant under test")
    ap.add_argument("--datasets", default=DEFAULT_DATASETS)
    ap.add_argument("--threads", default="16,64")
    ap.add_argument("--reps", type=int, default=3)
    ap.add_argument("--max-load", type=float, default=2.0)
    ap.add_argument("--settle", type=float, default=3.0)
    ap.add_argument("--timeout", type=float, default=2 * 3600)
    a = ap.parse_args()
    datasets, var = bench.manifest(), bench.variants()
    for v in (a.a, a.b):
        if v not in var["variants"]:
            sys.exit("unknown variant %s" % v)
    out_dir = os.path.join(bench.WORK, "ab")
    os.makedirs(out_dir, exist_ok=True)
    path = os.path.join(out_dir, "%s_vs_%s.jsonl" % (a.a, a.b))
    shas = {v: var["checkouts"][var["variants"][v]["checkout"]]["sha"] for v in (a.a, a.b)}
    done = {}
    if os.path.exists(path):
        for line in open(path):
            r = json.loads(line)
            if r["returncode"] == 0 and r["impl_sha"] == shas[r["variant"]]:
                done[r["key"]] = r
    run_opts = types.SimpleNamespace(systemd_scope=False, interval=0.25, timeout=a.timeout)
    jobs = []
    for d in a.datasets.split(","):
        if d not in datasets:
            print("skip %s: dataset not fetched" % d)
            continue
        for t in (int(x) for x in a.threads.split(",")):
            if t > bench.ncpu():
                print("skip %s t%d: only %d CPUs" % (d, t, bench.ncpu()))
                continue
            for rep in range(a.reps):
                for v in ((a.a, a.b) if rep % 2 == 0 else (a.b, a.a)):
                    j = dict(dataset=d, config="default", variant=v, threads=t, hashseed=0, rep=rep, suites=["ab"])
                    j["key"] = "ab|" + bench.job_key(j)
                    jobs.append(j)
    todo = [j for j in jobs if j["key"] not in done]
    print("%d A/B runs, %d to do" % (len(jobs), len(todo)))
    for i, j in enumerate(todo, 1):
        bench.wait_idle(a.max_load)
        time.sleep(a.settle)
        rec = bench.run_job(j, datasets, var, run_opts)
        with open(path, "a") as fh:
            fh.write(json.dumps(rec) + "\n")
        done[j["key"]] = rec
        print("[%d/%d] %-60s rc=%d %7.1f s" % (i, len(todo), j["key"], rec["returncode"], rec["wall_s"]), flush=True)
    summarize(jobs, done, a)


def summarize(jobs, done, a):
    groups = {}
    for j in jobs:
        r = done.get(j["key"])
        if r and r["returncode"] == 0:
            groups.setdefault((j["dataset"], j["threads"]), {}).setdefault(j["variant"], []).append(r)
    print("\n%-28s %4s %9s %9s %7s  %s" % ("dataset", "thr", a.a, a.b, "b/a", "outputs"))
    for (d, t), by in sorted(groups.items()):
        ra, rb = by.get(a.a, []), by.get(a.b, [])
        if not ra or not rb:
            continue
        ma, mb = (statistics.median(r["wall_s"] for r in x) for x in (ra, rb))
        hashes = {json.dumps(r["hashes"], sort_keys=True) for r in ra + rb}
        print("%-28s %4d %8.1fs %8.1fs %7.3f  %s" % (d, t, ma, mb, mb / ma, "same" if len(hashes) == 1 else "DIFFER"))
        phases = {}
        for v, rs in ((a.a, ra), (a.b, rb)):
            for r in rs:
                for p in r["phases"]:
                    phases.setdefault(p["phase"], {}).setdefault(v, []).append(p["dur_s"])
        for name, pv in phases.items():
            if a.a in pv and a.b in pv:
                print("      %-48s %8.2fs %8.2fs" % (name[:48], statistics.median(pv[a.a]), statistics.median(pv[a.b])))


if __name__ == "__main__":
    main()
