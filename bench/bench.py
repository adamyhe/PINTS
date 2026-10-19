#!/usr/bin/env python
"""PINTS benchmark runner: equivalence, runtime, peak memory and thread scaling.

    bench/run plan   --tier {smoke,standard,full} [--suites equivalence,timing,scaling] [--max-threads N]
    bench/run exec   [--suites ...] [--parallel N] [--systemd-scope] [--dry-run] [--retry-failed]
    bench/run status
    bench/run report [--out DIR]
    bench/run sysinfo

Every run executes the unmodified `scripts/pints_caller` of one implementation checkout (bench/variants.json)
in a fresh output directory and records wall time, CPU time, memory (summed PSS / RSS of the process tree,
system memory in use above the pre-run baseline, cgroup memory.peak when --systemd-scope is used),
`[phase]` markers (PINTS_PHASE_TIMING=1) and SHA-256 of the final peak files. Results are appended to
<work>/results/runs.jsonl; `exec` skips runs that already succeeded, so it can be interrupted and resumed.
"""
import argparse
import glob
import hashlib
import json
import os
import platform
import shutil
import socket
import subprocess
import sys
import threading
import time

HERE = os.path.dirname(os.path.abspath(__file__))
WORK = os.environ.get("PINTS_BENCH_WORK", os.path.join(os.path.dirname(HERE), "bench_work"))
RESULTS = os.path.join(WORK, "results")
RUNS_JSONL = os.path.join(RESULTS, "runs.jsonl")
PLAN = os.path.join(WORK, "plan.json")

CONFIGS = {"default": [], "sensitive": ["--sensitive"], "disable_small": ["--disable-small"]}
OUTPUT_GLOBS = ("*_bidirectional_peaks.bed", "*_divergent_peaks.bed", "*_unidirectional_peaks.bed")

T_SETS = ["T1_bw", "T1_bam", "T1_bam_ct", "T2_bam_k562", "T3_bw_gw_a549", "T4_bam_reps_hct116",
          "T5_bam_ct_k562_grocap"]
H_SETS = ["H1_bam_reps_karpas422", "H2_bw_caco2", "H3_bam_ct_gm12878_grocap", "H4_bw_liver", "H5_bw_nochr_caco2r2",
          "H6_bam_reps_bw_a549r2_karpas", "H7_batch_spleen", "H7_batch_testis", "H7_batch_kidney", "H7_batch_liver2",
          "H7_batch_cerebellum", "H7_batch_testis_b"]
ALL_SETS = T_SETS + H_SETS

# Tier definitions. Thread lists are capped at --max-threads (default: number of CPUs) when the plan is built.
TIERS = {
    "smoke": dict(
        eq_datasets=["T1_bw", "T1_bam"], eq_configs=["default"], eq_threads=8, determinism=[],
        timing_datasets=["T1_bw", "T1_bam"], timing_threads=[1, 8], reps=1,
        scaling_datasets=["T1_bam"], scaling_threads=[1, 2, 4, 8], scaling_reps=1, ref_scaling_threads=[1, 8]),
    "standard": dict(
        eq_datasets=ALL_SETS, eq_configs=list(CONFIGS), eq_threads=32,
        determinism=T_SETS,
        timing_datasets=["T2_bam_k562", "T3_bw_gw_a549", "T4_bam_reps_hct116", "T5_bam_ct_k562_grocap",
                         "H1_bam_reps_karpas422", "H2_bw_caco2", "H3_bam_ct_gm12878_grocap",
                         "H6_bam_reps_bw_a549r2_karpas"],
        timing_threads=[1, 16, 64], reps=3,
        scaling_datasets=["T3_bw_gw_a549", "T4_bam_reps_hct116", "H1_bam_reps_karpas422", "H3_bam_ct_gm12878_grocap"],
        scaling_threads=[1, 2, 4, 8, 16, 32, 64, "max"], scaling_reps=2, ref_scaling_threads=[1, 4, 16, 64]),
    "full": dict(
        eq_datasets=ALL_SETS, eq_configs=list(CONFIGS), eq_threads=32,
        determinism=T_SETS + ["H1_bam_reps_karpas422", "H2_bw_caco2", "H3_bam_ct_gm12878_grocap"],
        timing_datasets=ALL_SETS, timing_threads=[1, 16, 64, "max"], reps=3,
        scaling_datasets=["T2_bam_k562", "T3_bw_gw_a549", "T4_bam_reps_hct116", "T5_bam_ct_k562_grocap",
                          "H1_bam_reps_karpas422", "H2_bw_caco2", "H3_bam_ct_gm12878_grocap",
                          "H6_bam_reps_bw_a549r2_karpas"],
        scaling_threads=[1, 2, 4, 8, 16, 24, 32, 48, 64, "max"], scaling_reps=3,
        ref_scaling_threads=[1, 2, 4, 8, 16, 32, 64]),
}
EQ_VARIANTS = ["ref", "gate1", "gate2_off", "trackb", "gate2"]   # gate2 (fast math) is judged by drift, not hashes
# PR 4 (main + the empty-frame dtype fix) runs at reference speed, so it gets a reduced plan: equivalence on every
# dataset in the default config plus all configs on the dataset that triggers the fix, and timing on that dataset
PR4_VARIANT, PR4_ALL_CONFIGS, PR4_TIMING = "pr4", ["H3_bam_ct_gm12878_grocap"], ["H3_bam_ct_gm12878_grocap"]
REF_LIKE = ("ref", PR4_VARIANT)
TIMING_VARIANTS = ["ref", "gate1", "gate2"]
SCALING_VARIANTS = ["gate1", "gate2"]

# Rough single-run estimates (seconds) used only for `plan`: (ref 1 thread, ref 16, opt 1, opt 16) from the
# Biomni run on a 16-core machine; holdout values are guesses by data type. Expect +-2x.
EST = {
    "T1_bw": (80, 81, 18, 5.4), "T1_bam": (230, 229, 17, 5.8), "T1_bam_ct": (322, 322, 20, 6.9),
    "T2_bam_k562": (327, 286, 40, 10.6), "T3_bw_gw_a549": (1739, 443, 334, 44),
    "T4_bam_reps_hct116": (1930, 834, 244, 39), "T5_bam_ct_k562_grocap": (643, 455, 88, 19),
    "_gw_bw": (2400, 600, 450, 60), "_gw_bam": (6000, 2500, 800, 100), "_gw_bam_ct": (7000, 3000, 1200, 176),
}
# (ref, optimised) runtime relative to the default config; T2 (local) and T3 (Biomni) measurements
CONFIG_FACTOR = {"default": (1.0, 1.0), "sensitive": (1.05, 2.0), "disable_small": (1.9, 1.3)}


# ======================================================================================== helpers
def load_json(p, default=None):
    if os.path.exists(p):
        with open(p) as fh:
            return json.load(fh)
    return default


def variants():
    return load_json(os.path.join(HERE, "variants.json"))


def manifest():
    m = load_json(os.path.join(WORK, "manifest.json"))
    if not m:
        raise SystemExit("no %s; run `bench/run fetch` first" % os.path.join(WORK, "manifest.json"))
    return m["datasets"]


def ncpu():
    try:
        return len(os.sched_getaffinity(0))
    except AttributeError:
        return os.cpu_count() or 1


def job_key(j):
    return "%s|%s|%s|t%d|s%d|r%d%s" % (j["dataset"], j["config"], j["variant"], j["threads"], j["hashseed"], j["rep"],
                                       "|keep" if j.get("keep") else "")


def est_seconds(j):
    d = j["dataset"]
    if d in EST:
        e = EST[d]
    elif "_bw" in d:
        e = EST["_gw_bw"]
    elif "_ct_" in d:
        e = EST["_gw_bam_ct"]
    else:
        e = EST["_gw_bam"]
    ref1, ref16, opt1, opt16 = e
    t = j["threads"]
    cf_ref, cf_opt = CONFIG_FACTOR[j["config"]]
    if j["variant"] in REF_LIKE:
        if t >= 16:
            s = ref16
        else:  # reference scales poorly; interpolate in log2(threads)
            import math
            s = ref1 + (ref16 - ref1) * math.log2(t) / 4
        return s * cf_ref
    f = max(0.0, (opt16 / opt1 - 1 / 16) / (1 - 1 / 16))   # Amdahl serial fraction from the two measurements
    s = opt1 * (f + (1 - f) / t)
    if j["variant"] == "gate2":
        s *= 0.85
    return s * cf_opt


# ======================================================================================== plan
def build_plan(tier, suites, max_threads):
    T = TIERS[tier]

    def cap(ts):
        out = []
        for t in ts:
            t = max_threads if t == "max" else t
            if t <= max_threads and t not in out:
                out.append(t)
        return out

    jobs = {}

    def add(suite, exclusive, **j):
        j.setdefault("config", "default")
        j.setdefault("hashseed", 0)
        j.setdefault("rep", 0)
        k = job_key(j)
        if k in jobs:
            jobs[k]["suites"] = sorted(set(jobs[k]["suites"]) | {suite})
            jobs[k]["exclusive"] = jobs[k]["exclusive"] or exclusive
        else:
            j.update(key=k, suites=[suite], exclusive=exclusive)
            jobs[k] = j

    if "equivalence" in suites:
        te = min(T["eq_threads"], max_threads)
        for d in T["eq_datasets"]:
            for c in T["eq_configs"]:
                for v in EQ_VARIANTS:
                    add("equivalence", False, dataset=d, config=c, variant=v, threads=te,
                        keep=v in ("gate1", "gate2"))
        for d in T["determinism"]:   # thread count and hash seed must not change any output byte
            for t, seed in ((1, 0), (te, 1234), (max_threads, 99)):
                add("equivalence", False, dataset=d, variant="gate1", threads=t, hashseed=seed)
        if PR4_VARIANT in variants()["variants"]:
            for d in T["eq_datasets"]:
                for c in (list(CONFIGS) if d in PR4_ALL_CONFIGS else ["default"]):
                    add("equivalence", False, dataset=d, config=c, variant=PR4_VARIANT, threads=te)
    if "timing" in suites:
        if PR4_VARIANT in variants()["variants"]:
            for d in PR4_TIMING:
                if d not in T["timing_datasets"]:
                    continue
                for t in cap(T["timing_threads"]):
                    for rep in range(1 if t == 1 else T["reps"]):
                        add("timing", True, dataset=d, variant=PR4_VARIANT, threads=t, rep=rep)
        for d in T["timing_datasets"]:
            for t in cap(T["timing_threads"]):
                for rep in range(T["reps"]):
                    vs = TIMING_VARIANTS[rep % 3:] + TIMING_VARIANTS[:rep % 3]   # rotate order between repetitions
                    for v in vs:
                        if v == "ref" and t == 1 and rep > 0:
                            continue     # 1-thread reference runs take up to hours: one repetition
                        add("timing", True, dataset=d, variant=v, threads=t, rep=rep)
    if "scaling" in suites:
        for d in T["scaling_datasets"]:
            for t in cap(T["scaling_threads"]):
                for rep in range(T["scaling_reps"]):
                    for v in SCALING_VARIANTS:
                        add("scaling", True, dataset=d, variant=v, threads=t, rep=rep)
            for t in cap(T["ref_scaling_threads"]):
                add("scaling", True, dataset=d, variant="ref", threads=t)
    order = {"equivalence": 0, "timing": 1, "scaling": 2}
    out = sorted(jobs.values(), key=lambda j: (j["exclusive"], min(order[s] for s in j["suites"])))
    # keep the insertion order (alternation) inside each group: sorted() is stable
    measured = {}   # wall time of the last successful run of the same job, any commit (pr4 borrows ref's)
    for r in list(load_records(previous=True).values()) + list(load_records().values()):
        if r["returncode"] == 0:
            measured[r["key"]] = r["wall_s"]
    for j in out:
        alt = j["key"].replace("|%s|" % PR4_VARIANT, "|ref|") if j["variant"] == PR4_VARIANT else None
        m = measured.get(j["key"], measured.get(alt))
        j["est_s"] = round(m if m is not None else est_seconds(j), 1)
    return out


def cmd_plan(a):
    suites = a.suites.split(",")
    mt = a.max_threads or ncpu()
    jobs = build_plan(a.tier, suites, mt)
    if a.datasets:
        keep = set(a.datasets.split(","))
        jobs = [j for j in jobs if j["dataset"] in keep]
    os.makedirs(WORK, exist_ok=True)
    plan = {"tier": a.tier, "suites": suites, "max_threads": mt, "created": time.strftime("%Y-%m-%d %H:%M:%S"),
            "jobs": jobs}
    with open(PLAN, "w") as fh:
        json.dump(plan, fh, indent=1)
    done = completed_keys()
    print("plan: tier=%s  max_threads=%d  -> %s" % (a.tier, mt, PLAN))
    print("%-12s %6s %6s %10s" % ("suite", "jobs", "done", "est. hours"))
    for s in suites:
        js = [j for j in jobs if j["suites"][0] == s]
        print("%-12s %6d %6d %10.1f" % (s, len(js), sum(j["key"] in done for j in js),
                                        sum(j["est_s"] for j in js if j["key"] not in done) / 3600))
    tot = sum(j["est_s"] for j in jobs if j["key"] not in done) / 3600
    eq = sum(j["est_s"] for j in jobs if not j["exclusive"] and j["key"] not in done) / 3600
    print("total remaining: ~%.1f h sequential (equivalence part %.1f h; it can run with --parallel N)" % (tot, eq))
    print("Estimates use the last measured run of each job where there is one, otherwise figures scaled from a 16-core machine (+-2x).")


# ======================================================================================== sampling
def _read_kb(path, fields):
    out = {}
    try:
        with open(path) as fh:
            for line in fh:
                k = line.split(":", 1)[0]
                if k in fields:
                    out[k] = int(line.split()[1])
    except OSError:
        pass
    return out


def sys_used_bytes():
    if not os.path.exists("/proc/meminfo"):
        import psutil
        vm = psutil.virtual_memory()
        return vm.total - vm.available
    v = _read_kb("/proc/meminfo", ("MemTotal", "MemAvailable"))
    return (v["MemTotal"] - v["MemAvailable"]) * 1024


class TreeSampler(threading.Thread):
    """Samples the process tree rooted at `pid` every `interval` seconds.

    pss: summed proportional set size (shared pages split between the processes mapping them, so forked
    workers are not double counted); rss: summed RSS (overcounts shared pages); sys: MemTotal - MemAvailable
    minus the value before the run; cg: cgroup memory.current / memory.peak when the run has its own cgroup."""

    def __init__(self, pid, t0, interval):
        super().__init__(daemon=True)
        import psutil
        self.root = psutil.Process(pid)
        self.t0 = t0
        self.interval = interval
        self.base_sys = sys_used_bytes()
        self.samples = []   # (t, pss, rss, sys_delta, n_procs, cg_current)
        self.cg_dir = None
        self.cg_peak = None
        self._halt = threading.Event()
        self.linux = os.path.exists("/proc/self/smaps_rollup")

    def _cgroup(self):
        try:
            with open("/proc/%d/cgroup" % self.root.pid) as fh:
                for line in fh:
                    if line.startswith("0::"):
                        d = "/sys/fs/cgroup" + line.strip()[3:]
                        if os.path.exists(os.path.join(d, "memory.current")) and "pints-bench" in d:
                            return d
        except OSError:
            pass
        return None

    def run(self):
        import psutil
        while not self._halt.is_set():
            try:
                procs = [self.root] + self.root.children(recursive=True)
            except psutil.NoSuchProcess:
                break
            pss = rss = 0
            for p in procs:
                if self.linux:
                    v = _read_kb("/proc/%d/smaps_rollup" % p.pid, ("Pss", "Rss"))
                    pss += v.get("Pss", 0) * 1024
                    rss += v.get("Rss", 0) * 1024
                else:
                    try:
                        r = p.memory_info().rss
                    except (psutil.NoSuchProcess, psutil.AccessDenied):
                        continue
                    pss += r
                    rss += r
            if self.cg_dir is None:
                self.cg_dir = self._cgroup() or ""
            cg = None
            if self.cg_dir:
                try:
                    cg = int(open(os.path.join(self.cg_dir, "memory.current")).read())
                    self.cg_peak = int(open(os.path.join(self.cg_dir, "memory.peak")).read())
                except (OSError, ValueError):
                    pass
            self.samples.append((round(time.perf_counter() - self.t0, 3), pss, rss, sys_used_bytes() - self.base_sys,
                                 len(procs), cg))
            self._halt.wait(self.interval)

    def stop(self):
        self._halt.set()
        self.join()


# ======================================================================================== one run
def build_args(ds, config, outdir, threads):
    args = ["--save-to", outdir, "--file-prefix", "bench", "--thread", str(threads), "--dont-check-updates"]
    if ds.get("bam"):
        args += ["--bam-file"] + ds["bam"] + ["--exp-type", ds.get("exp_type", "PROcap")]
    else:
        args += ["--bw-pl"] + ds["bw_pl"] + ["--bw-mn"] + ds["bw_mn"]
    if ds.get("ct_bam"):
        args += ["--ct-bam"] + ds["ct_bam"]
    args += ds.get("extra_args", []) + CONFIGS[config]
    return args


def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def run_job(j, datasets, var, a):
    v = var["variants"][j["variant"]]
    impl = os.path.join(WORK, "impl", v["checkout"])
    outdir = os.path.join(WORK, "runs", j["key"].replace("|", "__"))
    if os.path.exists(outdir):
        shutil.rmtree(outdir)
    os.makedirs(outdir)
    ds = datasets[j["dataset"]]
    cmd = [sys.executable, os.path.join(HERE, "_pints_shim.py"), impl, "1" if j.get("keep") else "0"] + \
        build_args(ds, j["config"], outdir, j["threads"])
    if a.systemd_scope:
        cmd = ["systemd-run", "--user", "--scope", "--quiet", "-p", "MemoryAccounting=yes",
               "--unit=pints-bench-%d-%d" % (os.getpid(), int(time.time() * 1000) % 10 ** 9)] + cmd
    env = dict(os.environ)
    env.update({"PYTHONHASHSEED": str(j["hashseed"]), "OMP_NUM_THREADS": "1", "MKL_NUM_THREADS": "1",
                "OPENBLAS_NUM_THREADS": "1", "NUMBA_NUM_THREADS": "1", "LC_ALL": "C", "LANG": "C",
                "PYTHONWARNINGS": "ignore", "PINTS_PHASE_TIMING": "1", "TMPDIR": outdir})
    env.pop("PINTS_FAST_MATH", None)
    env.update(v["env"])
    load_before = os.getloadavg()[0]
    started = time.strftime("%Y-%m-%dT%H:%M:%S")
    t0 = time.perf_counter()
    so = open(os.path.join(outdir, "_stdout.log"), "w")
    se = open(os.path.join(outdir, "_stderr.log"), "w")
    proc = subprocess.Popen(cmd, env=env, cwd=outdir, stdout=so, stderr=subprocess.PIPE, text=True, bufsize=1)
    sampler = TreeSampler(proc.pid, t0, a.interval)
    sampler.start()
    markers = []

    def pump():
        for line in proc.stderr:
            se.write(line)
            if line.startswith("[phase]"):
                markers.append((round(time.perf_counter() - t0, 3), line[len("[phase]"):].rsplit(None, 1)[0].strip()))

    pt = threading.Thread(target=pump, daemon=True)
    pt.start()
    timed_out = False
    try:
        deadline = t0 + a.timeout
        while True:
            pid, status, ru = os.wait4(proc.pid, os.WNOHANG)
            if pid:
                break
            if time.perf_counter() > deadline:
                timed_out = True
                proc.kill()
            time.sleep(0.05)
    finally:
        wall = time.perf_counter() - t0
        sampler.stop()
        pt.join(timeout=10)
        so.close()
        se.close()
    rc = os.waitstatus_to_exitcode(status)
    proc.returncode = rc
    s = sampler.samples
    gb = 1e9
    phases, prev = [], 0.0
    for t, label in markers + [(wall, "tail (after last marker)")]:
        win = [x for x in s if prev < x[0] <= t]
        phases.append({"phase": label, "start": prev, "end": t, "dur_s": round(t - prev, 3),
                       "peak_pss_gb": max(x[1] for x in win) / gb if win else None,
                       "peak_sys_gb": max(x[3] for x in win) / gb if win else None})
        prev = t
    hashes = {}
    for g in OUTPUT_GLOBS:
        for p in sorted(glob.glob(os.path.join(outdir, g))):
            if "_sig_" not in os.path.basename(p):   # per-strand intermediates (only present in keep runs)
                hashes[os.path.basename(p)] = sha256(p)
    with open(os.path.join(outdir, "_samples.tsv"), "w") as fh:
        fh.write("t\tpss\trss\tsys_delta\tnprocs\tcgroup_current\n")
        for x in s:
            fh.write("\t".join("" if y is None else str(y) for y in x) + "\n")
    # large intermediates: per-chromosome coverage arrays; candidate tables are kept only for keep=True runs
    for p in glob.glob(os.path.join(outdir, "**", "*.npy"), recursive=True):
        os.remove(p)
    for p in glob.glob(os.path.join(outdir, "pybdt_*")) + glob.glob(os.path.join(outdir, "tmp*")):
        shutil.rmtree(p, ignore_errors=True)
    rec = {k: j[k] for k in ("key", "suites", "dataset", "config", "variant", "threads", "hashseed", "rep")}
    rec.update({
        "keep": bool(j.get("keep")), "impl_sha": var["checkouts"][v["checkout"]]["sha"], "host": socket.gethostname(),
        "started": started, "returncode": rc, "timed_out": timed_out, "wall_s": round(wall, 3),
        "cpu_s": round(ru.ru_utime + ru.ru_stime, 2), "max_single_rss_gb": ru.ru_maxrss * (1024 if sys.platform != "darwin" else 1) / gb,
        "peak_pss_gb": max((x[1] for x in s), default=0) / gb, "peak_rss_sum_gb": max((x[2] for x in s), default=0) / gb,
        "peak_sys_delta_gb": max((x[3] for x in s), default=0) / gb,
        "cgroup_peak_gb": sampler.cg_peak / gb if sampler.cg_peak else None,
        "max_procs": max((x[4] for x in s), default=0), "n_samples": len(s), "load_before": round(load_before, 2),
        "phases": phases, "hashes": hashes, "outdir": outdir, "cmd": " ".join(cmd),
    })
    rec["cpu_util"] = round(rec["cpu_s"] / wall, 2) if wall else None
    with open(os.path.join(outdir, "_record.json"), "w") as fh:
        json.dump(rec, fh, indent=1)
    if rc != 0:
        tail = open(os.path.join(outdir, "_stderr.log")).read()[-2000:]
        rec["stderr_tail"] = tail
    return rec


# ======================================================================================== exec / status
_lock = threading.Lock()


def append_record(rec):
    os.makedirs(RESULTS, exist_ok=True)
    with _lock, open(RUNS_JSONL, "a") as fh:
        fh.write(json.dumps(rec) + "\n")


def current_shas():
    var = variants()
    return {v: var["checkouts"][spec["checkout"]]["sha"] for v, spec in var["variants"].items()}


def load_records(previous=False):
    """Latest record per job key, counting only runs of the commit each variant is pinned to now in variants.json.

    Runs of other commits stay in runs.jsonl; with previous=True this returns the latest successful one per key
    instead (used by the report to compare a variant's new commit with the one it replaced)."""
    cur = current_shas()
    recs = {}
    if os.path.exists(RUNS_JSONL):
        with open(RUNS_JSONL) as fh:
            for line in fh:
                if line.strip():
                    r = json.loads(line)
                    if r["variant"] not in cur:
                        continue
                    if previous != (r.get("impl_sha") != cur[r["variant"]]):
                        continue
                    if previous and r["returncode"] != 0:
                        continue
                    recs[r["key"]] = r   # latest record per key wins
    return recs


def completed_keys():
    return {k for k, r in load_records().items() if r["returncode"] == 0}


def wait_idle(max_load, max_wait=600):
    t0 = time.time()
    while os.getloadavg()[0] > max_load and time.time() - t0 < max_wait:
        time.sleep(15)
    return os.getloadavg()[0]


def fmt_line(rec, i, n):
    return "[%d/%d] %-62s rc=%d %8.1fs pss=%6.2fGB sys=%6.2fGB cpu=%5.1fx" % (
        i, n, rec["key"], rec["returncode"], rec["wall_s"], rec["peak_pss_gb"], rec["peak_sys_delta_gb"],
        rec["cpu_util"] or 0)


def cmd_exec(a):
    plan = load_json(PLAN)
    if not plan:
        raise SystemExit("no plan; run `bench/run plan --tier ...` first")
    datasets, var = manifest(), variants()
    recs = load_records()
    jobs = plan["jobs"]
    if a.suites:
        want = set(a.suites.split(","))
        jobs = [j for j in jobs if want & set(j["suites"])]
    todo = []
    for j in jobs:
        r = recs.get(j["key"])
        if r and (r["returncode"] == 0 or not a.retry_failed):
            continue
        if j["dataset"] not in datasets:
            print("skip %s: dataset not fetched" % j["key"])
            continue
        if j["threads"] > ncpu():
            print("skip %s: %d threads > %d CPUs" % (j["key"], j["threads"], ncpu()))
            continue
        todo.append(j)
    print("%d runs to do (~%.1f h estimated sequential)" % (len(todo), sum(j["est_s"] for j in todo) / 3600))
    if a.dry_run:
        for j in todo:
            print("  %-70s ~%6.0fs" % (j["key"], j["est_s"]))
        return
    cmd_sysinfo(a, quiet=True)
    shared = [j for j in todo if not j["exclusive"]]
    excl = [j for j in todo if j["exclusive"]]
    n, i = len(todo), [0]
    failures = []
    streak = [0]

    class Abort(Exception):
        pass

    def one(j):
        if a.max_consecutive_failures and streak[0] >= a.max_consecutive_failures:
            raise Abort()
        rec = run_job(j, datasets, var, a)
        append_record(rec)
        with _lock:
            i[0] += 1
            print(fmt_line(rec, i[0], n), flush=True)
            if rec["returncode"] != 0:
                failures.append(rec["key"])
                streak[0] += 1
                print("    stderr: " + rec.get("stderr_tail", "")[-600:].replace("\n", "\n    "), flush=True)
            else:
                streak[0] = 0

    try:
        _exec_all(a, shared, excl, one)
    except Abort:
        print("stopping: %d consecutive failures (environment problem?). Fix it, then rerun with --retry-failed."
              % streak[0])
    print("done: %d runs, %d failed%s" % (i[0], len(failures), (": " + ", ".join(failures)) if failures else ""))


def _exec_all(a, shared, excl, one):
    if shared:
        par = max(1, a.parallel)
        # runs too wide to share the machine (e.g. the all-CPU determinism checks) go one at a time afterwards
        fits = [j for j in shared if par * j["threads"] <= ncpu()] if par > 1 else []
        alone = [j for j in shared if j not in fits]
        if par > 1 and alone:
            print("note: %d runs use more than %d/%d CPUs; they run one at a time after the others"
                  % (len(alone), ncpu(), par))
        if fits:
            from concurrent.futures import ThreadPoolExecutor
            with ThreadPoolExecutor(par) as ex:
                list(ex.map(one, fits))
        for j in alone:
            one(j)
    for j in excl:   # timing / scaling: strictly one at a time on an idle machine
        load = wait_idle(a.max_load)
        if load > a.max_load:
            print("warning: load average %.1f > %.1f before %s" % (load, a.max_load, j["key"]))
        time.sleep(a.settle)
        one(j)


def cmd_status(a):
    plan = load_json(PLAN)
    recs = load_records()
    if not plan:
        raise SystemExit("no plan")
    print("%-12s %6s %6s %6s %12s" % ("suite", "jobs", "ok", "failed", "est. left h"))
    for s in plan["suites"]:
        js = [j for j in plan["jobs"] if j["suites"][0] == s]
        ok = [j for j in js if recs.get(j["key"], {}).get("returncode") == 0]
        bad = [j for j in js if j["key"] in recs and recs[j["key"]]["returncode"] != 0]
        left = sum(j["est_s"] for j in js if j not in ok) / 3600
        print("%-12s %6d %6d %6d %12.1f" % (s, len(js), len(ok), len(bad), left))
        for j in bad:
            print("   failed: %s (rc=%s)" % (j["key"], recs[j["key"]]["returncode"]))


# ======================================================================================== sysinfo
def _sh(cmd):
    try:
        return subprocess.run(cmd, shell=True, capture_output=True, text=True, timeout=60).stdout.strip()
    except Exception as e:  # noqa: BLE001
        return "error: %s" % e


def cmd_sysinfo(a, quiet=False):
    import numpy
    info = {
        "time": time.strftime("%Y-%m-%dT%H:%M:%S"), "host": socket.gethostname(), "platform": platform.platform(),
        "python": sys.version, "ncpu_affinity": ncpu(), "cpu_count": os.cpu_count(),
        "lscpu": _sh("lscpu"), "meminfo": _sh("head -3 /proc/meminfo"), "numa": _sh("numactl -H 2>/dev/null | head -3"),
        "cgroup_cpu_max": _sh("cat /sys/fs/cgroup/cpu.max 2>/dev/null"),
        "cgroup_mem_max": _sh("cat /sys/fs/cgroup/memory.max 2>/dev/null"),
        "thp": _sh("cat /sys/kernel/mm/transparent_hugepage/enabled 2>/dev/null"),
        "governor": _sh("cat /sys/devices/system/cpu/cpu0/cpufreq/scaling_governor 2>/dev/null"),
        "filesystem": _sh("df -hT %s | tail -1" % WORK),
        "tools": {t: _sh("%s --version 2>&1 | head -1" % t) for t in ("bedtools", "samtools", "tabix", "bgzip")},
        "packages": _sh("%s -m pip list 2>/dev/null || uv pip list -p %s" % (sys.executable, sys.executable)),
        "numpy_runtime": _sh("%s -c 'import numpy; numpy.show_runtime()'" % sys.executable),
        "numpy": numpy.__version__,
        "impl": {k: {"sha": v["sha"], "head": _sh("git -C %s rev-parse HEAD" % os.path.join(WORK, "impl", k))}
                 for k, v in variants()["checkouts"].items()},
        "kit_commit": _sh("git -C %s rev-parse HEAD" % HERE),
    }
    for k, v in info["impl"].items():
        if v["sha"] != v["head"]:
            print("WARNING: checkout %s is at %s, variants.json pins %s" % (k, v["head"], v["sha"]))
    os.makedirs(RESULTS, exist_ok=True)
    hist = load_json(os.path.join(RESULTS, "sysinfo.json"), [])
    hist.append(info)
    with open(os.path.join(RESULTS, "sysinfo.json"), "w") as fh:
        json.dump(hist, fh, indent=1)
    if not quiet:
        print(json.dumps({k: info[k] for k in ("host", "platform", "ncpu_affinity", "meminfo", "cgroup_cpu_max",
                                               "cgroup_mem_max", "governor", "tools", "impl")}, indent=1))
        print(info["lscpu"][:1500])


# ======================================================================================== main
def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sp = ap.add_subparsers(dest="cmd", required=True)
    p = sp.add_parser("plan")
    p.add_argument("--tier", choices=tuple(TIERS), default="standard")
    p.add_argument("--suites", default="equivalence,timing,scaling")
    p.add_argument("--max-threads", type=int, default=0, help="cap thread counts (default: CPUs available)")
    p.add_argument("--datasets", default="", help="restrict to these dataset ids (comma-separated)")
    e = sp.add_parser("exec")
    e.add_argument("--suites", default="", help="only jobs belonging to these suites")
    e.add_argument("--parallel", type=int, default=1, help="concurrent runs for the equivalence suite only")
    e.add_argument("--systemd-scope", action="store_true",
                   help="run each job in its own systemd --user scope and record cgroup memory.peak")
    e.add_argument("--interval", type=float, default=0.25, help="memory sampling interval (s)")
    e.add_argument("--timeout", type=float, default=8 * 3600, help="per-run timeout (s)")
    e.add_argument("--max-load", type=float, default=2.0, help="wait for 1-min load below this before timed runs")
    e.add_argument("--settle", type=float, default=3.0, help="pause before each timed run (s)")
    e.add_argument("--retry-failed", action="store_true")
    e.add_argument("--max-consecutive-failures", type=int, default=5, help="stop after this many failures in a row "
                   "(0 = never)")
    e.add_argument("--dry-run", action="store_true")
    sp.add_parser("status")
    r = sp.add_parser("report")
    r.add_argument("--out", default=os.path.join(WORK, "report"))
    sp.add_parser("sysinfo")
    a = ap.parse_args()
    venv_py = os.environ.get("PINTS_BENCH_PYTHON")
    if a.cmd == "exec" and not (venv_py and os.path.realpath(venv_py) == os.path.realpath(sys.executable)):
        # A run started from another interpreter leaves PINTS children behind when the sampler fails to import.
        sys.exit("start runs with bench/run exec, which uses the kit's venv (%s is not it)" % sys.executable)
    if a.cmd == "plan":
        cmd_plan(a)
    elif a.cmd == "exec":
        cmd_exec(a)
    elif a.cmd == "status":
        cmd_status(a)
    elif a.cmd == "sysinfo":
        cmd_sysinfo(a)
    elif a.cmd == "report":
        sys.path.insert(0, HERE)
        import report
        report.main(a.out, load_records(), load_json(PLAN), load_json(os.path.join(RESULTS, "sysinfo.json"), []),
                    previous=load_records(previous=True), shas=current_shas())


if __name__ == "__main__":
    main()
