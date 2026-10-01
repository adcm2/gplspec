#!/usr/bin/env python3
"""Exercise the fail-closed dependency patch helper on disposable fixtures."""

from __future__ import annotations

import os
from pathlib import Path
import shutil
import subprocess
import tempfile

ROOT = Path(__file__).resolve().parents[1]
MODULE = ROOT / "cmake" / "ApplyPinnedOdrPatch.cmake"
BASE = {
    "include/a.h": "void first() {}\n",
    "include/b.h": "void second() {}\n",
}
PATCH_TEXT = """--- a/include/a.h
+++ b/include/a.h
@@ -1 +1 @@
-void first() {}
+inline void first() {}
--- a/include/b.h
+++ b/include/b.h
@@ -1 +1 @@
-void second() {}
+inline void second() {}
"""
PATCHED = {
    "include/a.h": "inline void first() {}\n",
    "include/b.h": "inline void second() {}\n",
}


def write_tree(root: Path, files: dict[str, str]) -> None:
    for relpath, contents in files.items():
        path = root / relpath
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(contents)


def snapshot(root: Path) -> dict[str, bytes]:
    return {relpath: (root / relpath).read_bytes() for relpath in BASE}


def invoke(cmake: str, script: Path, source: Path, patch: Path, build: Path,
           env: dict[str, str] | None = None) -> subprocess.CompletedProcess[str]:
    runner = script.parent / "invoke.cmake"
    runner.write_text(
        f'include("{MODULE}")\n'
        f'gplspec_apply_odr_patch("{source}" "{patch}" "{build}")\n'
    )
    return subprocess.run(
        [cmake, "-P", str(runner)], text=True, capture_output=True, env=env
    )


def fixture(tmp: Path, label: str, files: dict[str, str]) -> tuple[Path, Path, Path]:
    build = tmp / label / "build"
    source = build / "_deps" / "fake-dependency"
    source.mkdir(parents=True)
    write_tree(source, files)
    patch = tmp / label / "odr.patch"
    patch.write_text(PATCH_TEXT)
    return build, source, patch


def require(condition: bool, message: str) -> None:
    if not condition:
        raise AssertionError(message)


def main() -> None:
    cmake = shutil.which("cmake")
    require(cmake is not None, "cmake executable not found")
    with tempfile.TemporaryDirectory(prefix="gplspec-patch-guard-") as tmpdir:
        tmp = Path(tmpdir)

        build, source, patch = fixture(tmp, "pristine", BASE)
        result = invoke(cmake, tmp, source, patch, build)
        require(result.returncode == 0, f"pristine patch failed: {result.stderr}")
        require(snapshot(source) == {k: v.encode() for k, v in PATCHED.items()},
                "pristine source did not receive the complete patch")
        applied_snapshot = snapshot(source)
        result = invoke(cmake, tmp, source, patch, build)
        require(result.returncode == 0, f"already-applied patch failed: {result.stderr}")
        require(snapshot(source) == applied_snapshot,
                "already-applied check modified source")

        partial_files = dict(BASE)
        partial_files["include/a.h"] = PATCHED["include/a.h"]
        build, source, patch = fixture(tmp, "partial", partial_files)
        before = snapshot(source)
        result = invoke(cmake, tmp, source, patch, build)
        require(result.returncode != 0, "partial patch state was accepted")
        require(snapshot(source) == before,
                "partial patch rejection modified source")

        unexpected_files = dict(BASE)
        unexpected_files["include/a.h"] = "void changed() {}\n"
        build, source, patch = fixture(tmp, "unexpected", unexpected_files)
        before = snapshot(source)
        result = invoke(cmake, tmp, source, patch, build)
        require(result.returncode != 0, "unexpected source state was accepted")
        require(snapshot(source) == before,
                "unexpected source rejection modified source")

        build, source, patch = fixture(tmp, "git-failure", BASE)
        before = snapshot(source)
        fake_bin = tmp / "fake-bin"
        fake_bin.mkdir()
        fake_git = fake_bin / "git"
        fake_git.write_text("#!/bin/sh\nexit 73\n")
        fake_git.chmod(0o755)
        env = dict(os.environ)
        env["PATH"] = f"{fake_bin}{os.pathsep}{env.get('PATH', os.defpath)}"
        result = invoke(cmake, tmp, source, patch, build, env)
        require(result.returncode != 0, "git tool failure was accepted")
        require(snapshot(source) == before,
                "git tool failure modified source")

        build, source, patch = fixture(tmp, "outside-build", BASE)
        outside_source = tmp / "outside" / "dependency"
        write_tree(outside_source, BASE)
        before = snapshot(outside_source)
        result = invoke(cmake, tmp, outside_source, patch, build)
        require(result.returncode != 0, "out-of-build source was accepted")
        require(snapshot(outside_source) == before,
                "out-of-build rejection modified source")

    print("PASS: pristine, already-patched, partial, unexpected, git-failure, and out-of-build cases")


if __name__ == "__main__":
    main()
