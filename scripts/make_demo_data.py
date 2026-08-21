#!/usr/bin/env python3
"""Generate the synthetic demo dataset into demo/ (run from the repo root)."""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from sanger_seek.devtools.demogen import generate_demo  # noqa: E402

if __name__ == "__main__":
    target = Path(sys.argv[1]) if len(sys.argv) > 1 else Path(__file__).resolve().parents[1] / "demo"
    out = generate_demo(target)
    print(f"demo dataset written to {out}")
    for p in sorted(out.iterdir()):
        print(f"  {p.name}")
