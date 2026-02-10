---
title: Krylov-aware stochastic trace estimation
description: We introduce an algorithm for approximating the trace of a matrix function using careful deflation.
bibliography: krylov_aware_trace.bib
link-citations: true
link-bibliography: true
---


This is a companion piece to the paper: [https://arxiv.org/abs/2205.01736](https://arxiv.org/abs/2205.01736)

Code to replicate the figures in this document can be found on [Github](https://github.com/tchen-research/krylov_aware_trace).

## Introduction

Familiarity with  [Lanczos-FA](./lanczos-fa.html) or other Krylov subspace methods will be useful.  
An understanding of [Hutch++](https://ram900.hosting.nyu.edu/hutchplusplus/) is important for understanding the broader context, but not strictly necessary for this introduction.

Let $\vec{B}$ be a $d\times d$ symmetric matrix. 
An important primative in a number of state of the art trace estimation algorithms like Hutch++ [@meyer_musco_musco_woodruff_21] is finding an orthonormal matrix $\vec{Q}$ so that
$$
\vec{B} \approx \vec{Q} (\vec{Q}^\T \vec{B} \vec{Q}) \vec{Q}^\T.
$$
That is, finding the approximate span of the dominant eigenspace of $\vec{B}$.

A standard technique is by sketching [@halko_martinsson_tropp_11]. 
In particular, we can use the following algorithm:

:::{prf:algorithm}  Randomized range finder
**Input:** $\vec{A}\in\R^{n\times d}$
- sample Gaussian $d\times b$ matrix $\vec{\Omega}$
- compute $\vec{A}\vec{\Omega}$
- compute orthonormal basis $\vec{Q}$ for $\vec{A}\vec{\Omega}$

**Output:** $\vec{Q}$
:::

It can be shown that the $\vec{Q}$ produced by this algorithm is competitive with the best possible rank $b$ $\vec{Q}$.


## Krylov aware approximation

Often, $\vec{B} = f(\vec{A})$ for some matrix function $f(\vec{A})$.
In this case, it is standard to approximate products with $f(\vec{A})$ using a Krylov subspace methods, for instance the (block) Lanczos method for matrix function approximation.
Such methods approximate terms like $f(\vec{A})\vec{Z}$ by constructing the Krylov subspace 
$$
\mathcal{K}_{q+1}(\vec{A},\vec{Z}) = \operatorname{span}\{\vec{Z},\vec{A}\vec{Z}, \ldots, \vec{A}^q \vec{Z} \}.
$$

If we are using such a method with the randomized range finder, the first step is to approximate $f(\vec{A})\vec{\Omega}$ using $\mathcal{K}_{q+1}(\vec{A},\vec{\Omega})$. 
Then, we obtain an orthonormal basis $\vec{Q}$ for this space. 
However, in doing so, we implcitly construct the space $\mathcal{K}_{q+1}(\vec{A},\vec{\Omega}$, so we can obtain a better $\vec{Q}$ by using $\vec{\bar{Q}}_{q+1}$, an orthonormal basis for $\mathcal{K}_{q+1}(\vec{A},\vec{\Omega})$.

The ostensible downside is that to compute $\vec{Q}^\T f(\vec{A}) \vec{Q}$ now requires $q+1$ matrix-vector products with $\vec{A}$.
However, this is not actually the case.
The key observation made in this paper is the following:
If $\vec{\bar{Q}}_{q+1}$ is a basis for $\mathcal{K}_{q+1}(\vec{A},\vec{\Omega})$, then 
$$
\mathcal{K}_{n}(\vec{A},\vec{\bar{Q}}_{q+1})
= \mathcal{K}_{q+n}(\vec{A},\vec{\Omega}).
$$

This means that quantities such as $\vec{\bar{Q}}_{q+1}^\T f(\vec{A}) \vec{\bar{Q}}_{q+1}$ can be approximated from $\mathcal{K}_{q+n}(\vec{A},\vec{\Omega})$.
This requires just $nr$ additional matrix vector products with $\vec{A}$, rather than $n(q+1)$.

While this idea is simple, our algorithm can often perform much better than if matrix-vector products with $f(\vec{A})$ are treated as a black box. 
This is particularly true for functions like $f(x) = \sqrt{x}$, where we may require $q$ to be quite large to obtain a good approximation to $f(\vec{A}) \vec{\Omega}$, but even the basis obtained with $q=1$ (which is basically like sketching $\vec{A}$) works well. 

A similar idea was explored in [@persson_kressner_22] where it is shown that sketching $\vec{A}$ instead of $f(\vec{A})$ is a good idea for operator monotone functions.
