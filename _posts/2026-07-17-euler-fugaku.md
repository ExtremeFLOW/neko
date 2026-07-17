---
layout: post
title: Neko's compressible solver runs across 1.77 million cores of Fugaku
author: neko
---

Neko's new compressible flow solver has been exercised on up to 1,769,472 CPU cores of the supercomputer Fugaku at RIKEN R-CCS in Kobe, Japan — 36,864 of the machine's 158,976 nodes, roughly 23% of the system, and this represents the largest core count on which the code has been run to date.

Two pieces of new infrastructure made the runs possible: a hybrid MPI+OpenMP version of the code, and a new gather-scatter backend rebuilt on non-blocking MPI neighbourhood collectives. These new features will be included in the upcoming v1.1 release.
