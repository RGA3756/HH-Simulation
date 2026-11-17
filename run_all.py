"""
Reproduce every output of the project from scratch:

    python run_all.py

1. delete old outputs, record the computer and package versions
2. run the tests (stops here if any test fails)
3. run every simulation (sweep, failure temperatures, validation, refinement)
4. draw the figures and the GIF
5. build the interactive explorer and check its JavaScript against Python
6. write RESULTS_SUMMARY.md and METHODS_FILL_IN.md
"""
import re
import shutil
import subprocess
import sys
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parent
RESULTS, FIGURES = ROOT / "results", ROOT / "figures"


class Tee:
    """Print to the screen and to results/run_log.txt."""
    def __init__(self, path):
        self.file = open(path, "w")

    def __call__(self, msg=""):
        print(msg, flush=True)
        self.file.write(msg + "\n")
        self.file.flush()


def clean():
    for folder in (RESULTS, FIGURES):
        shutil.rmtree(folder, ignore_errors=True)
        folder.mkdir()
    for f in (ROOT / "web" / "explorer.html", ROOT / "RESULTS_SUMMARY.md", ROOT / "METHODS_FILL_IN.md"):
        f.unlink(missing_ok=True)


def main():
    t0 = time.time()
    clean()
    log = Tee(RESULTS / "run_log.txt")
    from hh_energy import animate, export_web, figures, report, sweep

    log("== 1. Environment")
    for k, v in report.write_env_info().items():
        log(f"  {k}: {v}")

    log("== 2. Tests")
    junit = RESULTS / "pytest_junit.xml"
    proc = subprocess.run([sys.executable, "-m", "pytest", "-q", f"--junitxml={junit}"], cwd=ROOT,
                          capture_output=True, text=True)
    # Remove the computer's network name from the test report (it is personal, not a result)
    junit.write_text(re.sub(r'hostname="[^"]*"', 'hostname=""', junit.read_text()))
    for line in proc.stdout.strip().splitlines()[-3:]:
        log("  " + line)
    if proc.returncode != 0:
        log(proc.stdout)
        log("A test failed: stopping before the experiments. Nothing below was run.")
        sys.exit(1)

    log("== 3. Experiments")
    sweep.run_experiments(RESULTS, log)

    log("== 4. Figures and animation")
    figures.make_all(log)
    animate.make_gif(log=log)

    log("== 5. Interactive explorer")
    export_web.build_html(log)
    js = export_web.js_consistency(log)
    if not js["pass"].all():
        log("WARNING: the JavaScript model does not match Python within 0.5 % (see results/js_consistency.csv)")

    log("== 6. Documents")
    report.write_all(log)
    log(f"Done in {time.time() - t0:.0f} s.")


if __name__ == "__main__":
    main()
