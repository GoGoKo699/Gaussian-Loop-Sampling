"""Run the documented deterministic and statistical verification suite."""
from pathlib import Path
import argparse
import json
import sys

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from check_identities import run as check_identities
from check_sampling import run as check_sampling
from check_topology import run as check_topology


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--samples", type=int, default=12000)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    if args.samples < 1000:
        parser.error("Use at least 1000 samples for the statistical diagnostics")
    results = {"identities": check_identities(), "sampling": check_sampling(args.samples),
               "topology": check_topology()}
    text = json.dumps(results, indent=2, allow_nan=False) + "\n"
    if args.output:
        args.output.write_text(text)
    print(text, end="")


if __name__ == "__main__":
    main()
