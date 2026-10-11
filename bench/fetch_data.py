#!/usr/bin/env python
"""Download, verify and prepare the benchmark datasets listed in bench/datasets.json.

    bench/run fetch [--only T2_bam_k562,H2_bw_caco2] [--parallel 4]

* ENCODE files are downloaded to <work>/data/raw and checked against the ENCODE md5.
* BAMs with `subset_chroms` are subset with samtools (full header kept, as PINTS allocates per @SQ).
* `strip_chr` bigWigs are rewritten with Ensembl-style chromosome names (H5).
* T1 datasets use the test data shipped in the reference checkout (tests/datasets).

Writes <work>/manifest.json (local paths, SHA-256 of every file fed to PINTS, and for BAMs a SHA-256 of
the alignment records alone). Files are compared with the SHA-256 recorded by the Biomni run; a mismatch
of a subset BAM is expected (samtools adds an @PG line with local paths, and compression differs between
builds) and is reported but not fatal. Use records_sha256 to compare inputs across machines.
"""
import argparse
import hashlib
import json
import os
import subprocess
import sys
import time
import urllib.request
from concurrent.futures import ThreadPoolExecutor

HERE = os.path.dirname(os.path.abspath(__file__))
WORK = os.environ.get("PINTS_BENCH_WORK", os.path.join(os.path.dirname(HERE), "bench_work"))
ENC = "https://www.encodeproject.org"
ROLES = ("bam", "ct_bam", "bw_pl", "bw_mn")


def _hash(path, algo):
    h = hashlib.new(algo)
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 22), b""):
            h.update(chunk)
    return h.hexdigest()


def fetch(acc, md5, raw):
    ext = None
    for e in (".bam", ".bigWig"):
        if os.path.exists(os.path.join(raw, acc + e)):
            ext = e
    if ext is None:
        req = urllib.request.Request("%s/files/%s/?format=json" % (ENC, acc), headers={"Accept": "application/json"})
        with urllib.request.urlopen(req, timeout=60) as r:
            meta = json.load(r)
        assert meta["md5sum"] == md5, "%s: ENCODE md5 changed (%s != %s)" % (acc, meta["md5sum"], md5)
        ext = {"bam": ".bam", "bigWig": ".bigWig"}[meta["file_format"]]
        dst = os.path.join(raw, acc + ext)
        tmp = dst + ".part"
        t0 = time.time()
        subprocess.check_call(["curl", "-fsSL", "--retry", "5", "--retry-delay", "10", "-C", "-", "-o", tmp,
                               ENC + meta["href"]])
        os.replace(tmp, dst)
        print("  downloaded %s (%.0f MB, %.0fs)" % (acc, os.path.getsize(dst) / 1e6, time.time() - t0), flush=True)
    dst = os.path.join(raw, acc + ext)
    got = _hash(dst, "md5")
    if got != md5:
        os.remove(dst)
        raise SystemExit("md5 mismatch for %s (got %s, expected %s); removed, rerun to retry" % (acc, got, md5))
    return dst


def samtools_index(path):
    if not os.path.exists(path + ".bai"):
        subprocess.check_call(["samtools", "index", "-@", "8", path])


def subset_bam(src, chroms, dst):
    if not os.path.exists(dst):
        samtools_index(src)
        subprocess.check_call(["samtools", "view", "--no-PG", "-@", "8", "-b", "-o", dst + ".part", src] + chroms)
        os.replace(dst + ".part", dst)
    samtools_index(dst)
    return dst


def bam_records_sha256(path):
    """SHA-256 of the alignment records as SAM text: independent of BAM compression and header @PG lines,
    so prepared inputs can be compared across machines and samtools builds."""
    h = hashlib.sha256()
    p = subprocess.Popen(["samtools", "view", "-@", "4", path], stdout=subprocess.PIPE)
    for chunk in iter(lambda: p.stdout.read(1 << 22), b""):
        h.update(chunk)
    if p.wait() != 0:
        raise RuntimeError("samtools view failed on %s" % path)
    return h.hexdigest()


def strip_chr_bw(src, dst):
    import pyBigWig
    if os.path.exists(dst):
        return dst
    bi = pyBigWig.open(src)
    chroms = [(c[3:] if c.startswith("chr") else c, L) for c, L in bi.chroms().items()]
    bo = pyBigWig.open(dst + ".part.bw", "w")
    bo.addHeader(chroms)
    for (c, _), (c2, _) in zip(bi.chroms().items(), chroms):
        iv = bi.intervals(c)
        if iv:
            bo.addEntries([c2] * len(iv), [x[0] for x in iv], ends=[x[1] for x in iv], values=[x[2] for x in iv])
    bo.close()
    bi.close()
    os.replace(dst + ".part.bw", dst)
    return dst


def prepare(ds, raw, proc, ref_checkout):
    ent = {k: ds[k] for k in ("id", "split", "description", "exp_type", "extra_args", "batch") if k in ds}
    ent["sha256"] = {}
    for role in ROLES:
        if "local" in ds and role in ds["local"]:
            ent[role] = [os.path.join(ref_checkout, p) for p in ds["local"][role]]
        elif role in ds.get("files", {}):
            paths = []
            for f in ds["files"][role]:
                p = fetch(f["accession"], f["md5"], raw)
                if role.endswith("bam"):
                    if ds.get("subset_chroms"):
                        p = subset_bam(p, ds["subset_chroms"], os.path.join(proc, "%s_%s.bam" % (ds["id"], f["accession"])))
                    else:
                        samtools_index(p)
                elif ds.get("strip_chr"):
                    p = strip_chr_bw(p, os.path.join(proc, "%s_%s.bw" % (ds["id"], f["accession"])))
                paths.append(p)
            ent[role] = paths
        else:
            continue
        for p in ent[role]:
            ent["sha256"][os.path.basename(p)] = _hash(p, "sha256")
            if p.endswith(".bam"):
                ent.setdefault("records_sha256", {})[os.path.basename(p)] = bam_records_sha256(p)
    want = ds.get("sha256_biomni", {})
    # Biomni names raw files <acc>.bam / <acc>.bigWig and prepared ones <id>_<acc>.bam / .bw, as here
    ent["sha256_matches_biomni"] = {k: (want[k] == v) for k, v in ent["sha256"].items() if k in want}
    return ent


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--only", default="", help="comma-separated dataset ids")
    ap.add_argument("--parallel", type=int, default=4, help="datasets prepared concurrently")
    a = ap.parse_args()
    spec = json.load(open(os.path.join(HERE, "datasets.json")))["datasets"]
    if a.only:
        keep = set(a.only.split(","))
        unknown = keep - {d["id"] for d in spec}
        if unknown:
            raise SystemExit("unknown dataset(s): %s" % ", ".join(sorted(unknown)))
        spec = [d for d in spec if d["id"] in keep]
    raw, proc = os.path.join(WORK, "data", "raw"), os.path.join(WORK, "data", "proc")
    os.makedirs(raw, exist_ok=True)
    os.makedirs(proc, exist_ok=True)
    ref_checkout = os.path.join(WORK, "impl", "ref")
    mpath = os.path.join(WORK, "manifest.json")
    manifest = json.load(open(mpath)) if os.path.exists(mpath) else {"datasets": {}}

    def one(ds):
        ent = prepare(ds, raw, proc, ref_checkout)
        bad = [k for k, ok in ent["sha256_matches_biomni"].items() if not ok]
        print("%-30s ok%s" % (ds["id"], ("  (sha256 differs from Biomni for %s)" % ", ".join(bad)) if bad else ""),
              flush=True)
        return ent

    failed = []
    with ThreadPoolExecutor(max(1, a.parallel)) as ex:
        futs = {ex.submit(one, d): d["id"] for d in spec}
        for fut, did in futs.items():
            try:
                manifest["datasets"][did] = fut.result()
            except BaseException as e:  # keep going; report at the end
                failed.append((did, e))
                print("%-30s FAILED: %s" % (did, e), file=sys.stderr, flush=True)
            json.dump(manifest, open(mpath, "w"), indent=1)
    print("manifest: %s (%d datasets)" % (mpath, len(manifest["datasets"])))
    if failed:
        raise SystemExit("%d dataset(s) failed: %s" % (len(failed), ", ".join(d for d, _ in failed)))


if __name__ == "__main__":
    main()
