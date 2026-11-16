### `pints.yulab.org` is returning 502 Bad Gateway

For the record, since this is what prompted the PR:

- **What:** `https://pints.yulab.org/` and the annotation API (`https://pints.yulab.org/api/ea/k562`) both return
  `502 Bad Gateway`, served by `nginx/1.20.1`. The front-end proxy is up, but the application behind it isn't
  answering.
- **When:** first seen on 2026-11-15, during the checks on the stack. Every check since has returned the same,
  most recently 2026-11-16 09:29 UTC.
- **Effect on `main` (and the 1.2.1 release):**
  - every `pints_caller --epig-annotation <biosample>` run stops after pairing with `FileNotFoundError`
    (`get_epig_annotation` logs `HTTP Error 502: Bad Gateway` and returns `""`, which
    `annotate_tre_with_epig_info` then cannot open). A local bigBed path fails the same way, because it is also
    sent to the server first;
  - `tests/functional_test.py` fails in `setUp` (`test_peak_calling`, `test_boundary_extension`,
    `test_bw_bam_consistency`), on every version back to `81b2dd9`, because it passes `epig_annotation="k562"`.
- **With this PR:** a local path works without the server, and a failed lookup gives a warning and peak calls
  without the annotation, so both the runs and the functional tests complete while the server is down.

The server itself needs a look on your side: the proxy is answering but its backend is not.

🤖 Generated with [Claude Code](https://claude.com/claude-code)
