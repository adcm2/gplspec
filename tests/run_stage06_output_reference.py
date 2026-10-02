#!/usr/bin/env python3
import hashlib
import pathlib
import subprocess
import sys
import tempfile


def read_manifest(path):
    result = {}
    for line in pathlib.Path(path).read_text().splitlines():
        digest, relative = line.split("  ", 1)
        result[relative] = digest
    return result


def main():
    if len(sys.argv) != 3:
        raise SystemExit("usage: run_stage06_output_reference.py EXECUTABLE MANIFEST")
    executable, manifest_path = sys.argv[1:]
    expected = read_manifest(manifest_path)
    with tempfile.TemporaryDirectory(prefix="gplspec-stage06-output-") as temp:
        root = pathlib.Path(temp)
        subprocess.run([executable, str(root)], check=True)
        actual = {}
        for path in sorted(root.rglob("*")):
            if path.is_file():
                actual[path.relative_to(root).as_posix()] = hashlib.sha256(
                    path.read_bytes()
                ).hexdigest()
    if actual != expected:
        missing = sorted(expected.keys() - actual.keys())
        extra = sorted(actual.keys() - expected.keys())
        changed = sorted(
            key for key in expected.keys() & actual.keys()
            if expected[key] != actual[key]
        )
        print(f"missing={missing}; extra={extra}; changed={changed}", file=sys.stderr)
        return 1
    print(f"PASS: {len(actual)} output files match the original-source baseline byte-for-byte")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
