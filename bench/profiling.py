#!/usr/bin/env python
"""Profiling runs of one implementation on a few datasets (default: gate1 on H1 and H3, default config).

    bench/run profile [--modes trace,pyspy,cprofile] [--datasets H1_...,H3_...] [--variant gate1] [--force]
    bench/run profile --dry-run

Modes (each is a separate run, so the profilers do not distort each other's timings):
  trace     function spans in the parent and every pool worker (see _pints_profile_shim.py), per-process PSS/RSS
            and I/O counters every 0.25 s, and every file that appears under the run directory (size, lifetime).
            Runs at 1, 16 and 64 threads on H1 and 1 and 64 on H3.
  pyspy     py-spy sampled stacks of all processes (--subprocesses --idle, raw collapsed format), 64 threads.
  cprofile  cProfile of the parent at 1 thread, where all the work runs in the parent (exact call counts).

Each run's peak files are hashed and compared with the benchmark records of the same variant, commit, dataset
and config, so a profiling run that changed the output would be flagged. Results go to <work>/profile/; the run
directories are deleted afterwards. Finally everything is packed into <work>/profile_results_<host>_<date>.tar.gz.
"""
import argparse
import glob
import json
import os
import shutil
import socket
import subprocess
import sys
import tarfile
import threading
import time

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import bench  # noqa: E402

PROFILE = os.path.join(bench.WORK, "profile")
RUNS = os.path.join(bench.WORK, "profile_runs")
H1, H3 = "H1_bam_reps_karpas422", "H3_bam_ct_gm12878_grocap"
THREADS = {"trace": {H1: [1, 16, 64], H3: [1, 64]}, "pyspy": {H1: [64], H3: [64]}, "cprofile": {H1: [1], H3: [1]}}
DEFAULT_THREADS = {"trace": [1, 64], "pyspy": [64], "cprofile": [1]}
EST_S = {(H1, 1): 530, (H1, 16): 73, (H1, 64): 55, (H3, 1): 663, (H3, 16): 91, (H3, 64): 74}
MODE_FACTOR = {"trace": 1.1, "pyspy": 1.3, "cprofile": 2.0}   # rough overhead guesses


class Sampler(threading.Thread):
    """Per-process memory and I/O of the process tree, plus an inventory of files under the run directory."""

    def __init__(self, pid, t0, outdir, interval=0.25, scan_every=1.0):
        super().__init__(daemon=True)
        import psutil
        self.root = psutil.Process(pid)
        self.t0, self.outdir, self.interval, self.scan_every = t0, outdir, interval, scan_every
        self.rows = []      # t, pid, ppid, pss, rss, rchar, wchar, read_bytes, write_bytes
        self.files = {}     # path -> [first_seen, last_seen, max_size]
        self._halt = threading.Event()

    @staticmethod
    def _io(pid):
        out = {}
        try:
            with open("/proc/%d/io" % pid) as fh:
                for line in fh:
                    k, v = line.split(":")
                    out[k] = int(v)
        except (OSError, ValueError):
            pass
        return out

    def _scan(self, t):
        for d, _, fs in os.walk(self.outdir):
            for f in fs:
                p = os.path.join(d, f)
                try:
                    sz = os.stat(p).st_size
                except OSError:
                    continue
                e = self.files.setdefault(os.path.relpath(p, self.outdir), [t, t, 0])
                e[1], e[2] = t, max(e[2], sz)

    def run(self):
        import psutil
        last_scan = -1e9
        while not self._halt.is_set():
            t = round(time.perf_counter() - self.t0, 3)
            try:
                procs = [self.root] + self.root.children(recursive=True)
            except psutil.NoSuchProcess:
                break
            for p in procs:
                m = bench._read_kb("/proc/%d/smaps_rollup" % p.pid, ("Pss", "Rss"))
                io = self._io(p.pid)
                try:
                    ppid = p.ppid()
                except psutil.NoSuchProcess:
                    continue
                self.rows.append((t, p.pid, ppid, m.get("Pss", 0) * 1024, m.get("Rss", 0) * 1024, io.get("rchar", 0),
                                  io.get("wchar", 0), io.get("read_bytes", 0), io.get("write_bytes", 0)))
            if t - last_scan >= self.scan_every:
                self._scan(t)
                last_scan = t
            self._halt.wait(self.interval)

    def stop(self):
        self._halt.set()
        self.join()
        self._scan(round(time.perf_counter() - self.t0, 3))


def plan(a):
    jobs = []
    for mode in a.modes.split(","):
        for ds in a.datasets.split(","):
            for t in THREADS[mode].get(ds, DEFAULT_THREADS[mode]):
                t = min(t, bench.ncpu())
                name = "%s__%s__t%d__%s" % (ds, a.variant, t, mode)
                if any(j["name"] == name for j in jobs):
                    continue
                est = EST_S.get((ds, t), 600) * MODE_FACTOR[mode]
                jobs.append({"mode": mode, "dataset": ds, "threads": t, "name": name, "est_s": est})
    return jobs


def done(job, sha):
    m = bench.load_json(os.path.join(PROFILE, job["name"], "meta.json"))
    return bool(m) and m.get("returncode") == 0 and m.get("impl_sha") == sha


def reference_hashes(variant, sha, dataset):
    """Peak-file hashes from benchmark runs of the same variant, commit, dataset and default config."""
    for r in bench.load_records().values():
        if (r["variant"], r.get("impl_sha"), r["dataset"], r["config"], r["returncode"]) == \
                (variant, sha, dataset, "default", 0) and r.get("hashes"):
            return r["hashes"]
    return None


def run(job, a, datasets, var):
    v = var["variants"][a.variant]
    impl = os.path.join(bench.WORK, "impl", v["checkout"])
    sha = var["checkouts"][v["checkout"]]["sha"]
    res = os.path.join(PROFILE, job["name"])
    outdir = os.path.join(RUNS, job["name"])
    for d in (res, outdir):
        shutil.rmtree(d, ignore_errors=True)
        os.makedirs(d)
    ds = datasets[job["dataset"]]
    mode = job["mode"]
    shim = [sys.executable, os.path.join(HERE, "_pints_profile_shim.py"), impl, "1" if mode == "trace" else "0",
            os.path.join(res, "cprofile.pstats") if mode == "cprofile" else "-"]
    cmd = shim + bench.build_args(ds, "default", outdir, job["threads"])
    if mode == "pyspy":
        pyspy = shutil.which("py-spy", path=os.path.dirname(sys.executable) + os.pathsep + os.environ.get("PATH", ""))
        if not pyspy:
            raise SystemExit("py-spy not found; rerun bash bench/setup.sh")
        cmd = [pyspy, "record", "--subprocesses", "--idle", "--nonblocking", "--rate", "100", "--format", "raw",
               "-o", os.path.join(res, "pyspy.txt"), "--"] + cmd
    env = dict(os.environ)
    env.update({"PYTHONHASHSEED": "0", "OMP_NUM_THREADS": "1", "MKL_NUM_THREADS": "1", "OPENBLAS_NUM_THREADS": "1",
                "NUMBA_NUM_THREADS": "1", "LC_ALL": "C", "LANG": "C", "PYTHONWARNINGS": "ignore",
                "PINTS_PHASE_TIMING": "1", "TMPDIR": outdir})
    env.pop("PINTS_FAST_MATH", None)
    env.update(v["env"])
    started = time.strftime("%Y-%m-%dT%H:%M:%S")
    load_before = os.getloadavg()[0]
    t0 = time.perf_counter()
    with open(os.path.join(res, "stdout.log"), "w") as so, open(os.path.join(res, "stderr.log"), "w") as se:
        proc = subprocess.Popen(cmd, env=env, cwd=outdir, stdout=so, stderr=se)
        try:
            sampler = Sampler(proc.pid, t0, outdir)
            sampler.start()
        except Exception:
            proc.kill()
            raise
        try:
            rc = proc.wait(timeout=a.timeout)
        except subprocess.TimeoutExpired:
            proc.kill()
            rc = proc.wait()
        wall = time.perf_counter() - t0
        sampler.stop()

    with open(os.path.join(res, "procs.tsv"), "w") as fh:
        fh.write("t\tpid\tppid\tpss\trss\trchar\twchar\tread_bytes\twrite_bytes\n")
        for r in sampler.rows:
            fh.write("\t".join(map(str, r)) + "\n")
    with open(os.path.join(res, "files.tsv"), "w") as fh:
        fh.write("path\tfirst_seen_s\tlast_seen_s\tmax_bytes\n")
        for p, (f, l, s) in sorted(sampler.files.items(), key=lambda kv: kv[1][0]):
            fh.write("%s\t%.2f\t%.2f\t%d\n" % (p, f, l, s))
    spans = parse_spans(os.path.join(res, "stderr.log"), t0)
    with open(os.path.join(res, "spans.tsv"), "w") as fh:
        fh.write("pid\tname\tcontext\tstart_s\tend_s\tdur_s\n")
        for s in spans:
            fh.write("%d\t%s\t%s\t%.4f\t%.4f\t%.4f\n" % s)

    hashes = {}
    for g in bench.OUTPUT_GLOBS:
        for p in sorted(glob.glob(os.path.join(outdir, g))):
            if "_sig_" not in os.path.basename(p):
                hashes[os.path.basename(p)] = bench.sha256(p)
    ref = reference_hashes(a.variant, sha, job["dataset"])
    meta = {"name": job["name"], "mode": mode, "dataset": job["dataset"], "variant": a.variant, "impl_sha": sha,
            "threads": job["threads"], "host": socket.gethostname(), "started": started, "load_before": load_before,
            "returncode": rc, "wall_s": round(wall, 3), "cmd": " ".join(cmd), "hashes": hashes,
            "matches_benchmark": None if ref is None or rc != 0 else hashes == ref,
            "peak_tree_pss_gb": round(max(tree_totals(sampler.rows).values() or [0]) / 1e9, 2)}
    with open(os.path.join(res, "meta.json"), "w") as fh:
        json.dump(meta, fh, indent=1)
    shutil.rmtree(outdir, ignore_errors=True)
    return meta


def parse_spans(path, t0):
    opened, out = {}, []
    with open(path, errors="replace") as fh:
        for line in fh:
            if not line.startswith("[span] "):
                continue
            try:
                head, ctx = line[7:].rstrip("\n").split("\t", 1)
                pid, kind, t, name = head.split(" ", 3)
                pid, t = int(pid), float(t) - t0
            except ValueError:
                continue
            key = (pid, name, ctx)
            if kind == "B":
                opened.setdefault(key, []).append(t)
            elif opened.get(key):
                start = opened[key].pop()
                out.append((pid, name, ctx, start, t, t - start))
    return sorted(out, key=lambda s: s[3])


def tree_totals(rows):
    tot = {}
    for r in rows:
        tot[r[0]] = tot.get(r[0], 0) + r[3]
    return tot


def pack():
    stamp = time.strftime("%Y%m%d_%H%M")
    out = os.path.join(bench.WORK, "profile_results_%s_%s.tar.gz" % (socket.gethostname().split(".")[0], stamp))
    with tarfile.open(out, "w:gz") as tf:
        tf.add(PROFILE, arcname="profile")
        si = os.path.join(bench.RESULTS, "sysinfo.json")
        if os.path.exists(si):
            tf.add(si, arcname="profile/sysinfo.json")
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--modes", default="trace,pyspy,cprofile")
    ap.add_argument("--datasets", default="%s,%s" % (H1, H3))
    ap.add_argument("--variant", default="gate1")
    ap.add_argument("--force", action="store_true", help="rerun jobs that already succeeded")
    ap.add_argument("--dry-run", action="store_true")
    ap.add_argument("--timeout", type=float, default=4 * 3600)
    a = ap.parse_args()
    venv_py = os.environ.get("PINTS_BENCH_PYTHON")
    if not (venv_py and os.path.realpath(venv_py) == os.path.realpath(sys.executable)):
        sys.exit("start this with bench/run profile, which uses the kit's venv (%s is not it)" % sys.executable)
    var = bench.variants()
    sha = var["checkouts"][var["variants"][a.variant]["checkout"]]["sha"]
    jobs = [j for j in plan(a) if a.force or not done(j, sha)]
    print("%d profiling runs to do (~%.0f min estimated)" % (len(jobs), sum(j["est_s"] for j in jobs) / 60))
    for j in jobs:
        print("  %-60s ~%5.0f s" % (j["name"], j["est_s"]))
    if a.dry_run:
        return
    datasets = bench.manifest()
    bench.cmd_sysinfo(a, quiet=True)
    for i, j in enumerate(jobs, 1):
        print("[%d/%d] %s ..." % (i, len(jobs), j["name"]), end=" ", flush=True)
        m = run(j, a, datasets, var)
        check = {True: "output matches benchmark", False: "OUTPUT DIFFERS FROM BENCHMARK",
                 None: "no benchmark hashes to compare"}[m["matches_benchmark"]]
        print("rc=%d  %.1f s  peak PSS %.1f GB  %s" % (m["returncode"], m["wall_s"], m["peak_tree_pss_gb"], check))
    print("results: %s" % pack())


if __name__ == "__main__":
    main()
