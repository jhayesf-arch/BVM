"""Every link between the web app's own files carries one version stamp (scripts/stamp_version.py),
so a publish never mixes new and cached files in a visitor's browser."""
import importlib.util
from pathlib import Path

spec = importlib.util.spec_from_file_location("stamp_version", Path(__file__).parent.parent / "scripts" / "stamp_version.py")
sv = importlib.util.module_from_spec(spec)
spec.loader.exec_module(sv)


def test_app_links_share_one_version_stamp():
    assert sv.check() == []
