#!/usr/bin/env python3
import pathlib
import shutil
import subprocess
import sys


def run(command, *, cwd=None):
    print("$", " ".join(map(str, command)), flush=True)
    subprocess.run(command, cwd=cwd, check=True)


source = pathlib.Path(sys.argv[1]).resolve()
build = pathlib.Path(sys.argv[2]).resolve()
deps = pathlib.Path(sys.argv[3]).resolve()
eigen = pathlib.Path(sys.argv[4]).resolve()
marker = build / ".gplspec-external-consumer-build"
if build.exists():
    if not marker.is_file() or marker.read_text() != str(source):
        raise RuntimeError(f"refusing to remove unowned external-consumer path: {build}")
    shutil.rmtree(build)
(build / "_deps").mkdir(parents=True)
marker.write_text(str(source))

dependency_names = {
    "numericconcepts": "NUMERICCONCEPTS",
    "fftwpp": "FFTWPP",
    "gaussquad": "GAUSSQUAD",
    "interpolation": "INTERPOLATION",
    "planetarymodel": "PLANETARYMODEL",
    "gshtrans": "GSHTRANS",
    "tomographymodels": "TOMOGRAPHYMODELS",
}
config = [
    "cmake", "-S", str(source / "tests/external_consumer"), "-B", str(build),
    "-G", "Ninja", "-DCMAKE_BUILD_TYPE=Release",
    f"-DGPLSPEC_SOURCE_DIR={source}",
    f"-DGPLSPEC_EIGEN_INCLUDE_DIR={eigen}",
]
for path in deps.glob("*-src"):
    key = path.name[:-4].lower()
    if key not in dependency_names:
        continue
    destination = build / "_deps" / path.name
    shutil.copytree(path, destination, ignore=shutil.ignore_patterns(".git"))
    config.append(f"-DFETCHCONTENT_SOURCE_DIR_{dependency_names[key]}={destination}")
if len(list((build / "_deps").glob("*-src"))) != len(dependency_names):
    raise RuntimeError("one or more pinned dependency sources are missing")
run(config)
run(["cmake", "--build", str(build), "--target", "external_consumer", "-j2"])
run([str(build / "external_consumer")])
