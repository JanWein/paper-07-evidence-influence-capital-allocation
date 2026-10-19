"""Reproduce the supplied synthetic analyses from the repository root."""
from pathlib import Path
import subprocess
import sys

ROOT = Path(__file__).resolve().parent
COMMANDS = ['code/partial_audit.py', 'code/reproduce_R2.py', 'code/reproduce_R1.py']
for script in COMMANDS:
    print(f"Running {script}", flush=True)
    subprocess.run([sys.executable, str(ROOT / script)], cwd=ROOT, check=True)
print("All selected scientific reproduction scripts completed successfully.")
