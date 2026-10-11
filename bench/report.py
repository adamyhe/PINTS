"""Build the benchmark report from <work>/results/runs.jsonl (called by `bench/run report`).

Writes <out>/report.md, CSV tables, PNG figures, and <work>/bench_results_<host>_<date>.tar.gz containing
everything needed to analyse the run elsewhere (records, sysinfo, plan, report; no peak files or data).
"""
import csv
import glob
import gzip
import json
import os
import statistics
import tarfile
import time
from collections import defaultdict

import numpy as np

KINDS = ("bidirectional", "divergent", "unidirectional")
HASH_VARIANTS = ("gate1", "gate2_off", "trackb")
TIMED = ("timing", "scaling")
# fixed categorical slots (reference palette, light mode): identity follows the implementation, never its rank
COLORS = {"ref": "#2a78d6", "gate1": "#eb6834", "gate2": "#1baf7a"}
LABELS = {"ref": "reference (main)", "gate1": "Gate 1", "gate2": "Gate 2 fast-math"}
INK, INK2, GRID = "#0b0b0b", "#52514e", "#e4e3df"


def _ok(r):
    return r["returncode"] == 0


def _hkey(r):
    return tuple(sorted(r["hashes"].items()))


def _write_csv(path, rows, cols):
    with open(path, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols, extrasaction="ignore")
        w.writeheader()
        for r in rows:
            w.writerow(r)


def _md_table(rows, cols, heads=None):
    heads = heads or cols
    out = ["| " + " | ".join(heads) + " |", "|" + "|".join("---" for _ in cols) + "|"]
    for r in rows:
        out.append("| " + " | ".join("" if r.get(c) is None else str(r.get(c)) for c in cols) + " |")
    return "\n".join(out)


def _fmt_s(x):
    if x is None:
        return None
    return "%.0f s" % x if x >= 100 else "%.1f s" % x


# ============================================================== equivalence
def equivalence(recs):
    by_dc = defaultdict(list)
    for r in recs.values():
        by_dc[(r["dataset"], r["config"])].append(r)
    rows, problems = [], []
    for (d, c), rs in sorted(by_dc.items()):
        refs = [r for r in rs if r["variant"] == "ref" and _ok(r) and r["hashes"]]
        ref_sets = {_hkey(r) for r in refs}
        row = {"dataset": d, "config": c, "ref_runs": len(refs),
               "ref_consistent": ("yes" if len(ref_sets) == 1 else ("NO (%d distinct)" % len(ref_sets))) if refs else "-"}
        canon = refs[0]["hashes"] if refs else None
        for v in HASH_VARIANTS + ("gate2",):
            vr = [r for r in rs if r["variant"] == v]
            if not vr:
                row[v] = ""
                continue
            failed = [r for r in vr if not _ok(r)]
            good = [r for r in vr if _ok(r)]
            if canon is None:
                row[v] = "no ref" + (" / %d failed" % len(failed) if failed else "")
                continue
            ident = [r for r in good if r["hashes"] == canon]
            cell = "%d/%d identical" % (len(ident), len(good))
            if failed:
                cell += ", %d FAILED" % len(failed)
            row[v] = cell
            if v in HASH_VARIANTS:
                for r in good:
                    if r["hashes"] != canon:
                        diff = sorted(k for k in set(canon) | set(r["hashes"]) if canon.get(k) != r["hashes"].get(k))
                        problems.append({"key": r["key"], "differing_files": ", ".join(diff)})
                for r in failed:
                    problems.append({"key": r["key"], "differing_files": "run failed (rc=%s)" % r["returncode"]})
        rows.append(row)
    return rows, problems


# ============================================================== gate 2 drift
def _candidates(rundir):
    rows = {}
    for f in sorted(glob.glob(os.path.join(rundir, "bench*_pl.bed.gz")) + glob.glob(os.path.join(rundir, "bench*_mn.bed.gz"))):
        if "_sig_" in f:
            continue
        with gzip.open(f, "rt") as fh:
            for line in fh:
                it = line.rstrip("\n").split("\t")
                rows[(os.path.basename(f), it[0], int(it[1]), int(it[2]), it[5])] = (float(it[4]), float(it[7]))
    return rows


def _final(rundir, kind):
    out = set()
    for f in glob.glob(os.path.join(rundir, "*_%s_peaks.bed" % kind)):
        if "_sig_" in f:
            continue
        for line in open(f):
            it = line.split("\t")
            out.add((os.path.basename(f), it[0], int(it[1]), int(it[2])))
    return out


def _rel(a, b):
    return np.abs(a - b) / np.maximum(np.maximum(np.abs(a), np.abs(b)), 1e-300)


def compare_pair(exact_dir, fast_dir, cutoff=0.1):
    """Biomni's pre-registered Gate 2 budget: identical candidate coordinates; p and padj max relative
    difference <= 1e-6; significance flips only within 1e-4 relative of the cutoff and <= 0.01% of significant
    candidates; Jaccard of every final call set >= 0.999."""
    ce, cf = _candidates(exact_dir), _candidates(fast_dir)
    res = {"n_candidates": len(ce), "coords_identical": set(ce) == set(cf)}
    keys = sorted(set(ce) & set(cf))
    if keys:
        pe, pf = np.array([ce[k] for k in keys]), np.array([cf[k] for k in keys])
        rp, rq = _rel(pe[:, 1], pf[:, 1]), _rel(pe[:, 0], pf[:, 0])
        res["pval_max_rel"] = float(np.max(np.where((pe[:, 1] < 1e-300) & (pf[:, 1] < 1e-300), 0, rp)))
        res["padj_max_rel"] = float(np.max(np.where((pe[:, 0] < 1e-300) & (pf[:, 0] < 1e-300), 0, rq)))
        se_, sf_ = pe[:, 0] < cutoff, pf[:, 0] < cutoff
        flips = np.where(se_ != sf_)[0]
        n_sig = int(se_.sum())
        res.update(n_significant=n_sig, n_flips=int(len(flips)),
                   flips_near_cutoff=bool(all(abs(pe[i, 0] - cutoff) / cutoff <= 1e-4 for i in flips)),
                   flip_fraction=(len(flips) / n_sig) if n_sig else 0.)
    for kind in KINDS:
        a, b = _final(exact_dir, kind), _final(fast_dir, kind)
        u = len(a | b)
        res["jaccard_" + kind] = (len(a & b) / u) if u else 1.
    res["pass"] = bool(keys and res["coords_identical"] and res["pval_max_rel"] <= 1e-6 and res["padj_max_rel"] <= 1e-6
                       and res["flips_near_cutoff"] and res["flip_fraction"] <= 1e-4
                       and all(res["jaccard_" + k] >= 0.999 for k in KINDS))
    return res


def drift(recs):
    pairs = defaultdict(dict)
    for r in recs.values():
        if r.get("keep") and _ok(r) and r["variant"] in ("gate1", "gate2"):
            pairs[(r["dataset"], r["config"], r["threads"])][r["variant"]] = r
    rows = []
    for (d, c, t), p in sorted(pairs.items()):
        if len(p) < 2:
            continue
        res = compare_pair(p["gate1"]["outdir"], p["gate2"]["outdir"])
        res.update(dataset=d, config=c, threads=t, wall_exact=p["gate1"]["wall_s"], wall_fast=p["gate2"]["wall_s"])
        rows.append(res)
    return rows


# ============================================================== timing / memory / scaling
def timed_groups(recs, config="default"):
    g = defaultdict(list)
    for r in recs.values():
        if _ok(r) and r["config"] == config and set(r["suites"]) & set(TIMED):
            g[(r["dataset"], r["threads"], r["variant"])].append(r)
    return g


def _med(rs, k):
    xs = [r[k] for r in rs if r.get(k) is not None]
    return statistics.median(xs) if xs else None


def timing_rows(g):
    rows = []
    for (d, t, v) in sorted(g):
        rs = g[(d, t, v)]
        ref = g.get((d, t, "ref"))
        ref1 = g.get((d, 1, "ref"))
        w = _med(rs, "wall_s")
        rows.append({
            "dataset": d, "threads": t, "variant": v, "n": len(rs), "wall_median_s": round(w, 2),
            "wall_min_s": round(min(r["wall_s"] for r in rs), 2), "wall_max_s": round(max(r["wall_s"] for r in rs), 2),
            "speedup_vs_ref_same_threads": round(_med(ref, "wall_s") / w, 2) if ref and v != "ref" else None,
            "speedup_vs_ref_1_thread": round(_med(ref1, "wall_s") / w, 2) if ref1 else None,
            "cpu_util": round(_med(rs, "cpu_util"), 1) if _med(rs, "cpu_util") else None,
            "peak_pss_gb": round(_med(rs, "peak_pss_gb"), 2),
            "peak_sys_delta_gb": round(_med(rs, "peak_sys_delta_gb"), 2),
            "cgroup_peak_gb": round(_med(rs, "cgroup_peak_gb"), 2) if _med(rs, "cgroup_peak_gb") else None,
            "max_single_rss_gb": round(_med(rs, "max_single_rss_gb"), 2),
            "pss_ratio_vs_ref": round(_med(rs, "peak_pss_gb") / _med(ref, "peak_pss_gb"), 2)
            if ref and v != "ref" and _med(ref, "peak_pss_gb") else None,
        })
    return rows


def phase_rows(g):
    """Median phase durations for the optimised variants (PINTS_PHASE_TIMING markers)."""
    rows = []
    for (d, t, v), rs in sorted(g.items()):
        if v == "ref":
            continue
        acc = defaultdict(list)
        mem = defaultdict(list)
        for r in rs:
            dur, peak = defaultdict(float), {}
            for p in r.get("phases", []):   # a label can repeat (e.g. once per sample): sum within a run
                dur[p["phase"]] += p["dur_s"]
                if p.get("peak_pss_gb") is not None:
                    peak[p["phase"]] = max(peak.get(p["phase"], 0), p["peak_pss_gb"])
            for ph, x in dur.items():
                acc[ph].append(x)
                if ph in peak:
                    mem[ph].append(peak[ph])
        total = _med(rs, "wall_s")
        for ph, ds in acc.items():
            m = statistics.median(ds)
            rows.append({"dataset": d, "threads": t, "variant": v, "phase": ph, "dur_median_s": round(m, 2),
                         "share": round(m / total, 3) if total else None,
                         "peak_pss_gb": round(statistics.median(mem[ph]), 2) if mem[ph] else None})
    return rows


# ============================================================== figures
def _style(ax, title, xlabel, ylabel):
    ax.set_title(title, loc="left", fontsize=10, color=INK)
    ax.set_xlabel(xlabel, fontsize=9, color=INK2)
    ax.set_ylabel(ylabel, fontsize=9, color=INK2)
    ax.tick_params(colors=INK2, labelsize=8, length=0)
    ax.grid(True, color=GRID, linewidth=0.8)
    ax.set_axisbelow(True)
    for s in ax.spines.values():
        s.set_visible(False)


def _thread_axis(ax, threads):
    from matplotlib.ticker import FixedLocator, NullLocator
    ax.set_xscale("log", base=2)
    ax.xaxis.set_major_locator(FixedLocator(sorted(threads)))
    ax.xaxis.set_minor_locator(NullLocator())
    ax.set_xticklabels([str(t) for t in sorted(threads)])


def scaling_figures(g, figdir):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    out = []
    datasets = sorted({d for (d, t, v) in g if v == "gate1"})
    for d in datasets:
        series = {}
        for v in ("ref", "gate1", "gate2"):
            pts = sorted((t, _med(rs, "wall_s"), _med(rs, "peak_pss_gb")) for (dd, t, vv), rs in g.items() if dd == d and vv == v)
            if pts:
                series[v] = pts
        if len(series.get("gate1", [])) < 3:
            continue
        fig, axes = plt.subplots(1, 3, figsize=(13, 3.6), dpi=150)
        fig.patch.set_facecolor("#fcfcfb")
        for ax in axes:
            ax.set_facecolor("#fcfcfb")
        # 1) wall time, log-log, with ideal linear scaling from 1 thread
        ax = axes[0]
        for v, pts in series.items():
            ts, ws = [p[0] for p in pts], [p[1] for p in pts]
            ax.plot(ts, ws, color=COLORS[v], linewidth=2, marker="o", markersize=4, label=LABELS[v])
            if ts[0] == 1 and v != "ref":
                ax.plot(ts, [ws[0] / t for t in ts], color=COLORS[v], linewidth=1, linestyle=":", alpha=0.7)
        allt = sorted({p[0] for pts in series.values() for p in pts})
        _thread_axis(ax, allt)
        ax.set_yscale("log")
        _style(ax, "%s: wall time (dotted = ideal)" % d, "threads", "seconds (log)")
        ax.legend(fontsize=7, frameon=False)
        # 2) self-relative speedup and parallel efficiency share one quantity: speedup over own 1-thread run
        ax = axes[1]
        for v, pts in series.items():
            if pts[0][0] != 1:
                continue
            ts = [p[0] for p in pts]
            ax.plot(ts, [pts[0][1] / p[1] for p in pts], color=COLORS[v], linewidth=2, marker="o", markersize=4,
                    label=LABELS[v])
        tmax = max(p[0] for pts in series.values() for p in pts)
        ax.plot([1, tmax], [1, tmax], color=INK2, linewidth=1, linestyle=":", label="ideal")
        _thread_axis(ax, allt)
        ax.set_yscale("log", base=2)
        from matplotlib.ticker import ScalarFormatter
        ax.yaxis.set_major_formatter(ScalarFormatter())
        _style(ax, "speedup over own 1-thread run", "threads", "speedup (log2)")
        ax.legend(fontsize=7, frameon=False)
        # 3) peak memory (summed PSS of the process tree)
        ax = axes[2]
        for v, pts in series.items():
            ax.plot([p[0] for p in pts], [p[2] for p in pts], color=COLORS[v], linewidth=2, marker="o", markersize=4,
                    label=LABELS[v])
        _thread_axis(ax, allt)
        ax.set_ylim(bottom=0)
        _style(ax, "peak memory (summed PSS)", "threads", "GB")
        ax.legend(fontsize=7, frameon=False)
        fig.tight_layout()
        p = os.path.join(figdir, "scaling_%s.png" % d)
        fig.savefig(p, facecolor=fig.get_facecolor())
        plt.close(fig)
        out.append(p)
    return out


def speedup_figure(trows, figdir):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    threads = sorted({r["threads"] for r in trows if r["speedup_vs_ref_same_threads"]})
    if not threads:
        return []
    out = []
    for t in threads:
        rs = [r for r in trows if r["threads"] == t and r["speedup_vs_ref_same_threads"]]
        ds = sorted({r["dataset"] for r in rs})
        vs = [v for v in ("gate1", "gate2") if any(r["variant"] == v for r in rs)]
        fig, ax = plt.subplots(figsize=(max(6, 0.9 * len(ds) + 2), 3.4), dpi=150)
        fig.patch.set_facecolor("#fcfcfb")
        ax.set_facecolor("#fcfcfb")
        w = 0.8 / len(vs)
        for i, v in enumerate(vs):
            ys = [next((r["speedup_vs_ref_same_threads"] for r in rs if r["dataset"] == d and r["variant"] == v), 0)
                  for d in ds]
            xs = np.arange(len(ds)) + (i - (len(vs) - 1) / 2) * w
            ax.bar(xs, ys, width=w * 0.92, color=COLORS[v], label=LABELS[v])
            for x, y in zip(xs, ys):
                if y:
                    ax.annotate("%.1f×" % y, (x, y), xytext=(0, 2), textcoords="offset points", ha="center",
                                fontsize=6, color=INK2)
        ax.set_xticks(np.arange(len(ds)))
        ax.set_xticklabels([d if d.startswith("T1") else d.split("_")[0] + "_" + d.split("_")[-1] for d in ds],
                           rotation=30, ha="right")
        _style(ax, "Speedup over the reference at %d threads (median wall time)" % t, "", "× faster than reference")
        ax.legend(fontsize=7, frameon=False)
        fig.tight_layout()
        p = os.path.join(figdir, "speedup_t%d.png" % t)
        fig.savefig(p, facecolor=fig.get_facecolor())
        plt.close(fig)
        out.append(p)
    return out


# ============================================================== main
def main(out, recs, plan, sysinfo):
    for r in recs.values():   # older records of keep runs also hashed the per-strand _sig_ intermediates
        r["hashes"] = {k: v for k, v in r.get("hashes", {}).items() if "_sig_" not in k}
    os.makedirs(out, exist_ok=True)
    figdir = os.path.join(out, "figures")
    os.makedirs(figdir, exist_ok=True)
    md = ["# PINTS optimisation benchmark", "",
          "Generated %s from %d run records." % (time.strftime("%Y-%m-%d %H:%M"), len(recs)), ""]
    if sysinfo:
        si = sysinfo[-1]
        cpu = next((l.split(":", 1)[1].strip() for l in si.get("lscpu", "").splitlines() if l.startswith("Model name")), "?")
        md += ["## Machine", "",
               "- host `%s`, %s" % (si["host"], si["platform"]),
               "- CPU: %s; %s CPUs available to the run" % (cpu, si["ncpu_affinity"]),
               "- memory: %s" % si.get("meminfo", "").splitlines()[0] if si.get("meminfo") else "- memory: ?",
               "- cgroup cpu.max `%s`, memory.max `%s`, governor `%s`, THP `%s`" % (
                   si.get("cgroup_cpu_max") or "-", si.get("cgroup_mem_max") or "-", si.get("governor") or "-",
                   si.get("thp") or "-"),
               "- implementations: " + ", ".join("%s `%s`" % (k, v["sha"][:7]) for k, v in si["impl"].items()), ""]
    if plan:
        md += ["Plan: tier `%s`, suites %s, max threads %s." % (plan["tier"], ", ".join(plan["suites"]),
                                                              plan["max_threads"]), ""]
        missing = [j["key"] for j in plan["jobs"] if j["key"] not in recs]
        if missing:
            md += ["**%d planned runs have no record yet**; tables below are partial." % len(missing), ""]

    # equivalence
    erows, problems = equivalence(recs)
    _write_csv(os.path.join(out, "equivalence.csv"), erows,
               ["dataset", "config", "ref_runs", "ref_consistent"] + list(HASH_VARIANTS) + ["gate2"])
    n_cmp = sum(int(c.split("/")[1].split()[0]) for r in erows for v in HASH_VARIANTS
                for c in [r.get(v) or ""] if "identical" in c)
    md += ["## 1. Equivalence (byte identity of final peak files vs the reference)", "",
           "Every successful run of `gate1`, `gate2_off` (Gate 2 branch with the flag unset) and `trackb` "
           "(default flags) must reproduce the reference's bidirectional, divergent and unidirectional peak files "
           "byte for byte, at any thread count and hash seed. `gate2` (fast math) is shown for information "
           "only; it is judged by the drift budget in section 2.", "",
           "**%s** (%d variant runs compared)." % ("All identical" if not problems else "%d mismatches or failures"
                                                    % len(problems), n_cmp), "",
           _md_table(erows, ["dataset", "config", "ref_consistent"] + list(HASH_VARIANTS) + ["gate2"],
                     ["dataset", "config", "ref self-consistent", "gate1", "gate2 (flag off)", "trackB", "gate2 fast"]), ""]
    if problems:
        md += ["### Mismatches / failures", "", _md_table(problems, ["key", "differing_files"]), ""]

    # drift
    drows = drift(recs)
    dcols = ["dataset", "config", "threads", "n_candidates", "coords_identical", "pval_max_rel", "padj_max_rel",
             "n_flips", "flip_fraction", "jaccard_bidirectional", "jaccard_divergent", "jaccard_unidirectional",
             "wall_exact", "wall_fast", "pass"]
    _write_csv(os.path.join(out, "drift.csv"), drows, dcols)
    if drows:
        for r in drows:
            for k in ("pval_max_rel", "padj_max_rel"):
                r[k + "_fmt"] = "%.1e" % r[k] if r.get(k) is not None else None
            for k in KINDS:
                r["j_" + k[:3]] = "%.5f" % r["jaccard_" + k]
            r["gain"] = "%.1f%%" % (100 * (1 - r["wall_fast"] / r["wall_exact"]))
        md += ["## 2. Gate 2 drift (PINTS_FAST_MATH=1 vs Gate 1)", "",
               "Budget fixed before testing: identical candidates; max relative Δp and Δpadj ≤ 1e-6; flips only "
               "within 1e-4 of the cutoff and ≤ 0.01%% of significant candidates; Jaccard ≥ 0.999 for every final set. "
               "**%d/%d pass.** Wall-time gains here come from equivalence runs, which may have shared the machine; "
               "use section 3 for timing." % (sum(r["pass"] for r in drows), len(drows)), "",
               _md_table(drows, ["dataset", "config", "n_candidates", "coords_identical", "pval_max_rel_fmt",
                                 "padj_max_rel_fmt", "n_flips", "j_bid", "j_div", "j_uni", "gain", "pass"],
                         ["dataset", "config", "candidates", "coords same", "max rel Δp", "max rel Δpadj", "flips",
                          "J bidir", "J div", "J unidir", "wall gain", "pass"]), ""]

    # timing / memory
    g = timed_groups(recs)
    trows = timing_rows(g)
    tcols = ["dataset", "threads", "variant", "n", "wall_median_s", "wall_min_s", "wall_max_s",
             "speedup_vs_ref_same_threads", "speedup_vs_ref_1_thread", "cpu_util", "peak_pss_gb",
             "peak_sys_delta_gb", "cgroup_peak_gb", "max_single_rss_gb", "pss_ratio_vs_ref"]
    _write_csv(os.path.join(out, "timing_memory.csv"), trows, tcols)
    figs = speedup_figure(trows, figdir) + scaling_figures(g, figdir)
    if trows:
        for r in trows:
            r["wall"] = _fmt_s(r["wall_median_s"]) + (" (n=%d)" % r["n"] if r["n"] > 1 else "")
            r["spd"] = ("%.1f×" % r["speedup_vs_ref_same_threads"]) if r["speedup_vs_ref_same_threads"] else ""
            r["spd1"] = ("%.1f×" % r["speedup_vs_ref_1_thread"]) if r["speedup_vs_ref_1_thread"] and r["variant"] != "ref" else ""
            r["mem"] = "%.2f" % r["peak_pss_gb"]
            r["memr"] = ("%.2f×" % r["pss_ratio_vs_ref"]) if r["pss_ratio_vs_ref"] else ""
        md += ["## 3. Runtime and peak memory (default config)", "",
               "Timed runs ran one at a time, alternating implementations between repetitions. Speedup is the ratio "
               "of median wall times. Memory is the peak of the summed PSS over the process tree (shared pages "
               "split between the forked workers); `timing_memory.csv` also has summed RSS, system memory above "
               "baseline, cgroup peak and the largest single-process RSS.", ""]
        for p in figs:
            if "speedup_t" in p:
                md += ["![](figures/%s)" % os.path.basename(p), ""]
        md += [_md_table(trows, ["dataset", "threads", "variant", "wall", "spd", "spd1", "cpu_util", "mem", "memr"],
                         ["dataset", "threads", "variant", "wall (median)", "speedup vs ref", "vs ref 1 thread",
                          "CPU util (cores)", "peak PSS GB", "PSS vs ref"]), ""]
    sc = [p for p in figs if "scaling_" in p]
    if sc:
        md += ["## 4. Thread scaling", "",
               "Left: median wall time (dotted lines: ideal scaling from each implementation's own 1-thread run). "
               "Middle: speedup over the implementation's own 1-thread run. Right: peak summed PSS.", ""]
        md += ["![](figures/%s)" % os.path.basename(p) for p in sc]
        md += [""]
    prows = phase_rows(g)
    _write_csv(os.path.join(out, "phases.csv"), prows, ["dataset", "threads", "variant", "phase", "dur_median_s",
                                                        "share", "peak_pss_gb"])
    if prows:
        tmax = max(r["threads"] for r in prows)
        sel = [r for r in prows if r["variant"] == "gate1" and r["threads"] in (1, tmax)]
        sel.sort(key=lambda r: (r["dataset"], r["threads"], -r["dur_median_s"]))
        md += ["## 5. Where the time goes (Gate 1, phase markers)", "",
               "Median duration and peak PSS per phase at 1 and %d threads, from `PINTS_PHASE_TIMING` markers timestamped "
               "by the runner; the first phase also includes interpreter start-up and imports. Full table: "
               "`phases.csv`." % tmax, "",
               _md_table(sel, ["dataset", "threads", "phase", "dur_median_s", "share", "peak_pss_gb"],
                         ["dataset", "threads", "phase", "seconds", "share of wall", "peak PSS GB"]), ""]
    fails = [r for r in recs.values() if not _ok(r)]
    if fails:
        md += ["## Failed runs", "", _md_table([{"key": r["key"], "rc": r["returncode"],
                                                "tail": (r.get("stderr_tail") or "").strip().splitlines()[-1:] or [""]}
                                               for r in fails], ["key", "rc", "tail"]), ""]
    with open(os.path.join(out, "report.md"), "w") as fh:
        fh.write("\n".join(md))
    work = os.environ.get("PINTS_BENCH_WORK", os.path.dirname(out))
    host = sysinfo[-1]["host"] if sysinfo else "host"
    tgz = os.path.join(work, "bench_results_%s_%s.tar.gz" % (host, time.strftime("%Y%m%d_%H%M")))
    with tarfile.open(tgz, "w:gz") as tf:
        for p in ("results", "plan.json", "manifest.json"):
            if os.path.exists(os.path.join(work, p)):
                tf.add(os.path.join(work, p), arcname=os.path.join("bench_results", p))
        tf.add(out, arcname="bench_results/report")
        # per-run memory traces and phase logs are small and useful for follow-up analysis
        for d in glob.glob(os.path.join(work, "runs", "*")):
            for f in ("_samples.tsv", "_record.json", "_stderr.log"):
                p = os.path.join(d, f)
                if os.path.exists(p) and os.path.getsize(p) < 20e6:
                    tf.add(p, arcname=os.path.join("bench_results", "runs", os.path.basename(d), f))
    print("report: %s\nbundle: %s" % (os.path.join(out, "report.md"), tgz))
