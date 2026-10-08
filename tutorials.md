---
title: "Tutorials"
layout: default
permalink: /tutorials/
---

# Scientific tutorials

The tutorial collection documents model construction, geometry/mapping, tomography, and output interpretation. These pages are explanatory material; compiling or reading a tutorial is not a validation of its scientific examples.

## Available tutorials

{% for tutorial in site.tutorials %}
- [{{ tutorial.title }}]({{ tutorial.url | relative_url }}) — {{ tutorial.description }}
{% endfor %}

For executable prerequisites and current CMake usage, see [Getting started](getting-started.md). For the input, parameter, output, plot-consumer, and validation map, see [Scientific workflows](benchmark.md).
