#!/usr/bin/env python3
"""Extract the C++ example fence from a Markdown tutorial for compilation."""
import pathlib
import re
import sys

source, destination = map(pathlib.Path, sys.argv[1:])
text = source.read_text(encoding="utf-8")
match = re.search(r"```cpp\s*\n(.*?)\n```", text, re.DOTALL)
if not match:
    raise SystemExit(f"no C++ code fence found in {source}")
destination.parent.mkdir(parents=True, exist_ok=True)
destination.write_text(match.group(1) + "\n", encoding="utf-8")
