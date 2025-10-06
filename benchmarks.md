---
title: "Benchmarks"
layout: default
permalink: /benchmarks/
---


# Benchmarks

This section contains performance benchmarks and validation tests for gplspec.

## Available Benchmarks

{% for benchmark in site.benchmarks %}
- [{{ benchmark.title }}]({{ benchmark.url | relative_url }}) - {{ benchmark.description }}
{% endfor %}

## Running Benchmarks

Instructions on how to compile and run the benchmarks:

```bash
cd build
make clean_bench_1
./clean_bench_1
```