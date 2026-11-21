"""Run scan_runner.py jobs from a plan file with bounded concurrency; skips jobs whose summary exists.

usage: python results/scan_queue.py plan.txt --jobs 6
plan.txt: one scan_runner argument string per line ('#' comments). Logs: results/scan_logs/<root>__<tag>.out|err
"""
import argparse
import hashlib
import shlex
import subprocess
import sys
import time
from pathlib import Path

HERE = Path(__file__).resolve().parent
LOGS = HERE / "scan_logs"


def _name(args):
    toks = shlex.split(args)
    root = toks[toks.index("--root") + 1]
    name = "_".join(t.lstrip("-") for t in toks if t != "--root" and t != root)
    if len(name) > 160:  # keep log paths under the Windows MAX_PATH limit
        name = name[:150] + "_" + hashlib.sha1(name.encode()).hexdigest()[:8]
    return root, name


def _ok(out_path):
    """A job succeeded only if scan_runner printed an 'ok' RESULT line (exit code 0 is not enough)."""
    try:
        return any(ln.startswith("RESULT") and " ok " in ln for ln in Path(out_path).read_text().splitlines())
    except OSError:
        return False


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("plan")
    ap.add_argument("--jobs", type=int, default=6)
    a = ap.parse_args()
    lines = [ln.strip() for ln in Path(a.plan).read_text().splitlines() if ln.strip() and not ln.startswith("#")]
    LOGS.mkdir(exist_ok=True)
    pending, running, marked = list(lines), [], set()
    while pending or running:
        for p, done in running:
            if p.poll() is not None and done not in marked:
                marked.add(done)
                if p.returncode == 0 and _ok(done.with_suffix(".out")):
                    done.write_text("ok")
                else:
                    print(time.strftime("%H:%M:%S"), "FAILED", done.stem, flush=True)
        running = [(p, n) for p, n in running if p.poll() is None]
        while pending and len(running) < a.jobs:
            args = pending.pop(0)
            root, name = _name(args)
            done = LOGS / f"{root}__{name}.done"
            if done.exists() or _ok(done.with_suffix(".out")):
                continue
            out = open(LOGS / f"{root}__{name}.out", "w")
            err = open(LOGS / f"{root}__{name}.err", "w")
            cmd = [sys.executable, "-W", "ignore", str(HERE / "scan_runner.py")] + shlex.split(args)
            p = subprocess.Popen(cmd, stdout=out, stderr=err, cwd=HERE.parent)
            running.append((p, done))
            print(time.strftime("%H:%M:%S"), "start", root, name, flush=True)
        time.sleep(10)
    print("all done", flush=True)


if __name__ == "__main__":
    main()
