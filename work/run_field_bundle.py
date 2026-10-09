"""Run a producer and bind its field bundles to unchanged input bytes."""

import argparse
import hashlib
import json
from pathlib import Path
import subprocess


def sha256(path):
    digest = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--workdir", required=True)
    parser.add_argument("--input", action="append", default=[])
    parser.add_argument("--bundle", action="append", required=True)
    parser.add_argument("command", nargs=argparse.REMAINDER)
    args = parser.parse_args()
    command = args.command[1:] if args.command and args.command[0] == "--" else args.command
    if not command:
        parser.error("a producer command is required after --")
    workdir = Path(args.workdir).resolve()
    inputs = {name: sha256(workdir / name) for name in args.input}
    subprocess.run(command, cwd=workdir, check=True)
    if any(sha256(workdir / name) != digest for name, digest in inputs.items()):
        raise RuntimeError("a declared input changed while the producer was running")
    for bundle_name in args.bundle:
        path = workdir / bundle_name / "manifest.json"
        manifest = json.loads(path.read_text(encoding="utf-8"))
        for identity in manifest.get("inputs", []):
            name = identity.get("path")
            if name not in inputs:
                raise RuntimeError(
                    f"bundle {bundle_name} declares input not listed to runner: {name}"
                )
            identity["sha256"] = inputs[name]
        temporary = path.with_suffix(".json.tmp")
        temporary.write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n",
                             encoding="utf-8")
        temporary.replace(path)


if __name__ == "__main__":
    main()
