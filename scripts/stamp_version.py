"""Stamp a version query (?v=...) on every link between the web app's own files.

GitHub Pages lets browsers cache files for about 10 minutes, so after a publish a normal reload can
mix new HTML with old JavaScript. Changing the query string on each publish makes browsers fetch the
new files. Run this before committing a change to app/ (it is idempotent):

    python scripts/stamp_version.py            # version = UTC timestamp
    python scripts/stamp_version.py --check    # exit 1 if any link is unstamped or stamps differ
"""
from __future__ import annotations

import argparse
import re
import sys
from datetime import datetime, timezone
from pathlib import Path

APP = Path(__file__).resolve().parent.parent / "app"
FILES = sorted(APP.glob("*.js")) + sorted(APP.glob("*.html"))
# a quoted relative reference to one of the app's own .js/.css files, with or without a stamp
LINK = re.compile(r'(["\'])((?:\./)?[A-Za-z0-9_-]+\.(?:js|css))(?:\?v=[\w.-]*)?\1')


def _own(name: str) -> bool:
    return (APP / name.removeprefix("./")).is_file()


def stamp(version: str) -> int:
    n = 0
    for f in FILES:
        text = f.read_text()
        new = LINK.sub(lambda m: f"{m[1]}{m[2]}?v={version}{m[1]}" if _own(m[2]) else m[0], text)
        if new != text:
            f.write_text(new)
            n += 1
    return n


def check() -> list[str]:
    found: dict[str, set[str]] = {}
    for f in FILES:
        for m in LINK.finditer(f.read_text()):
            if _own(m[2]):
                v = re.search(r"\?v=([\w.-]*)", m[0])
                found.setdefault(v[1] if v else "<none>", set()).add(f.name)
    return [] if len(found) == 1 and "<none>" not in found else [f"{v}: {sorted(fs)}" for v, fs in found.items()]


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--check", action="store_true")
    ap.add_argument("--version")
    a = ap.parse_args()
    if a.check:
        bad = check()
        print("\n".join(bad) or "all links stamped with one version")
        sys.exit(1 if bad else 0)
    v = a.version or datetime.now(timezone.utc).strftime("%Y%m%d%H%M")
    print(f"stamped v={v} in {stamp(v)} file(s)")
