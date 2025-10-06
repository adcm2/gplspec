---
title: "Benchmarks"
layout: default
---

# Benchmarks

This section contains performance benchmarks and validation tests for gplspec.

## Available Benchmarks

{% for benchmarks in site.benchmarks %}
- [{{ benchmarks.title }}]({{ benchmarks.url | relative_url }}) - {{ benchmarks.description }}
{% endfor %}

## Running Benchmarks

Instructions on how to compile and run the benchmarks:

```bash
cd build
make clean_bench_1
./clean_bench_1
```