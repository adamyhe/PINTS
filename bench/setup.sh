#!/usr/bin/env bash
# One-time setup for the PINTS benchmark kit (Linux x86_64).
#
#   bash bench/setup.sh [WORK_DIR]        # default WORK_DIR: ./bench_work (or $PINTS_BENCH_WORK)
#
# Creates, under WORK_DIR:
#   venv/          uv virtualenv (CPython 3.12) with the pinned requirements
#   tools/bin/     bgzip, tabix, samtools (htslib/samtools 1.22.1, built from source) and bedtools 2.31.0,
#                  unless all four are already on PATH (set PINTS_BENCH_BUILD_TOOLS=1 to build anyway)
#   impl/<name>/   git worktrees of the implementations listed in bench/variants.json
#   env.sh         sourced by bench/run
set -euo pipefail

KIT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO="$(dirname "$KIT")"
WORK="$(mkdir -p "${1:-${PINTS_BENCH_WORK:-$REPO/bench_work}}" && cd "${1:-${PINTS_BENCH_WORK:-$REPO/bench_work}}" && pwd)"
HTS_VER=1.22.1
BEDTOOLS_VER=2.31.0   # newest release with a static Linux binary (2.31.1 ships source only)
JOBS="$(nproc 2>/dev/null || echo 4)"

say() { printf '\n==> %s\n' "$*"; }

command -v uv >/dev/null || { echo "uv not found; install it first: https://docs.astral.sh/uv/getting-started/installation/" >&2; exit 1; }
command -v git >/dev/null || { echo "git not found" >&2; exit 1; }

# ---------------------------------------------------------------- python
say "Python environment ($WORK/venv)"
[ -x "$WORK/venv/bin/python" ] || uv venv -q -p 3.12 "$WORK/venv"
uv pip install -q -p "$WORK/venv/bin/python" -r "$KIT/requirements.txt"
"$WORK/venv/bin/python" -c "import numpy, scipy, pandas, pysam, pyBigWig, pybedtools, numba, psutil; print('python deps ok, numpy', numpy.__version__, 'numba', numba.__version__)"

# ---------------------------------------------------------------- command-line tools
mkdir -p "$WORK/tools/bin" "$WORK/tools/src"
export PATH="$WORK/tools/bin:$PATH"
have_all=1
for t in bgzip tabix samtools bedtools; do command -v "$t" >/dev/null || have_all=0; done
if [ "$have_all" = 1 ] && [ "${PINTS_BENCH_BUILD_TOOLS:-0}" != 1 ]; then
  say "Using bgzip/tabix/samtools/bedtools already on PATH"
else
  for t in gcc make curl; do command -v "$t" >/dev/null || { echo "$t is required to build htslib/samtools" >&2; exit 1; }; done
  cd "$WORK/tools/src"
  for pkg in htslib samtools; do
    if [ ! -x "$WORK/tools/bin/$([ $pkg = htslib ] && echo bgzip || echo samtools)" ]; then
      say "Building $pkg $HTS_VER"
      curl -fsSL -o "$pkg-$HTS_VER.tar.bz2" "https://github.com/samtools/$pkg/releases/download/$HTS_VER/$pkg-$HTS_VER.tar.bz2"
      tar xjf "$pkg-$HTS_VER.tar.bz2"
      (cd "$pkg-$HTS_VER" &&
        if [ $pkg = htslib ]; then
          ./configure -q --prefix="$WORK/tools" --disable-libcurl --disable-gcs --disable-s3 --disable-bz2 --disable-lzma
        else
          ./configure -q --prefix="$WORK/tools" --with-htslib="$WORK/tools/src/htslib-$HTS_VER" --without-curses \
            --disable-bz2 --disable-lzma
        fi && make -s -j"$JOBS" && make -s install) >"$pkg.build.log" 2>&1 ||
        { echo "building $pkg failed (needs zlib headers, e.g. zlib1g-dev / zlib-devel); see $WORK/tools/src/$pkg.build.log" >&2; exit 1; }
    fi
  done
  if [ ! -x "$WORK/tools/bin/bedtools" ] && command -v bedtools >/dev/null; then
    ln -s "$(command -v bedtools)" "$WORK/tools/bin/bedtools"
  elif [ ! -x "$WORK/tools/bin/bedtools" ]; then
    [ "$(uname -s)-$(uname -m)" = Linux-x86_64 ] ||
      { echo "no static bedtools binary for $(uname -s)-$(uname -m); install bedtools and rerun" >&2; exit 1; }
    say "Downloading bedtools $BEDTOOLS_VER (static binary)"
    curl -fsSL -o "$WORK/tools/bin/bedtools" \
      "https://github.com/arq5x/bedtools2/releases/download/v$BEDTOOLS_VER/bedtools.static"
    chmod +x "$WORK/tools/bin/bedtools"
  fi
fi
for t in bgzip tabix samtools bedtools; do printf '  %-9s %s\n' "$t" "$(command -v $t)"; done
bedtools --version; samtools --version | head -1

# ---------------------------------------------------------------- implementations
say "Implementation checkouts ($WORK/impl)"
cd "$REPO"
[ "$(git rev-parse --is-shallow-repository)" = true ] && git fetch -q --unshallow origin
mkdir -p "$WORK/impl"
"$WORK/venv/bin/python" - "$KIT/variants.json" <<'EOF' | while read -r name ref sha; do
import json, sys
for k, v in json.load(open(sys.argv[1]))["checkouts"].items():
    print(k, v["ref"], v["sha"])
EOF
  if ! git cat-file -e "$sha^{commit}" 2>/dev/null; then
    git fetch -q origin "$ref" || git fetch -q https://github.com/adamyhe/PINTS.git "$ref"
  fi
  git cat-file -e "$sha^{commit}" || { echo "commit $sha ($ref) not available" >&2; exit 1; }
  dst="$WORK/impl/$name"
  if [ -d "$dst" ]; then
    [ "$(git -C "$dst" rev-parse HEAD)" = "$sha" ] || git -C "$dst" checkout -q --detach "$sha"
  else
    git worktree add -q --detach "$dst" "$sha"
  fi
  printf '  %-7s %s (%s)\n' "$name" "${sha:0:7}" "$ref"
done

cat >"$WORK/env.sh" <<EOF
export PINTS_BENCH_WORK="$WORK"
export PATH="$WORK/tools/bin:\$PATH"
export PINTS_BENCH_PYTHON="$WORK/venv/bin/python"
EOF

say "Done. Next:  bench/run fetch      (downloads ~5.8 GB of ENCODE data into $WORK/data)"
