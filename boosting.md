---
title: Median trick in high dimension
description: Short description of the median trick in high dimensions.
keywords: median trick, boosting, randomized algorithms, high dimensions
link-citations: true
---

## Introduction

The "median trick" (sometimes called "boosting") is a clever way of turning a (scaler valued) randomized algorithm that succeeds with some constant probability (e.g. 2/3) into one that succeeds with arbitrarily high probability $1-\delta$ with only logarithmic overhead.

The idea is simple: suppose we have a randomized algorithm that produces an approximation $x$ to $x^*$ satisfying $\mathbb{P}[ |x^*-x| < \varepsilon ] > 2/3$.
We can boost the success probability by running the algorithm multiple times and taking the median of the results.

1. Run the algorithm $m$ times to get estimates $x_1, x_2, \ldots, x_m$ 
2. Let $\widehat{x} = \operatorname{median}(x_1, \ldots, x_m)$

:::{prf:proposition}
 For $m = O(\log(1/\delta))$, we have that $\mathbb{P}[ |x^*-\widehat{x}| < \varepsilon ] > 1-\delta$.
:::

Notably, the dependence on the target failure probability $\delta$ is only *logarithmic*.

The proof is simple and standard, so we will just provide a high level overview.

:::{prf:proof}
:class: dropdown
:enumerated: false
By definition of the median, half of the estimates $x_i$ are below $\widehat{x}$ and half are above.
Thus, the only way for $\widehat{x}$ to be outside of $[x^*-\varepsilon, x^*+\varepsilon]$ is if at least half of the estimates $x_i$ are outside of this interval.
However, since each $x_i$ live within this interval with probability at least $2/3$, it is (exponentially) unlikely that half of them live outside of the interval.
:::

## High dimensions

What about in high dimensions?
I.e. if we have a randomized algorithm that produces an approximation $\vec{x}$ to $\widehat{\vec{x}}\in\mathbb{R}^d$ satisfying $\mathbb{P}[ \|\widehat{\vec{x}}-\vec{x}\| < \varepsilon ] > 2/3$?
While the median trick in 1d appears in almost all course notes for a class on randomized algorithms, finding a statement of a high-dimensional analog is surprisingly much [more difficult](https://cstheory.stackexchange.com/questions/32536/generalizing-the-median-trick-to-higher-dimensions). 


The simplest approach is to apply the one-dimensional estimator along each dimension.
The main issue with this is that $\|\vec{x} - \widehat{\vec{x}}\|$ being small is not necessarily the same thing as each coordinate being close (i.e. $\|\vec{x} - \widehat{\vec{x}}\|_\infty$ being small).
So we will lose something in the accuracy from converting between norms.
For instance, if $\| \cdot \|$ is the standard Euclidan norm, then we will lose a factor of $\sqrt{d}$ in the accuracy.


### A better estimator

Here is a better approach that only loses a factor of 3 on the accuracy. 
I was introduced to this approach by [Chris Musco](https://www.chrismusco.com/), but I'd expect its a folklore type result.


1.  Make independent copies $\vec{x}_1, \ldots, \vec{x}_m$ of $\vec{x}$
1. For all $i,j\in[m]$, define $d(i,j) = \| \vec{x}_i - \vec{x}_j \|$
1. For all $i\in[m]$, define $c(i) = \operatorname{median}(d(i,1),\ldots, d(i,t))$
1. $i^* = \operatorname{argmin}_{i\in[m]} c(i)$
1. $\widehat{\vec{x}} = \vec{x}_{i^*}$

:::{prf:proposition}
For $m = O(\log(1/\delta))$, we have that $\mathbb{P}[ |\vec{x}^*-\widehat{\vec{x}}| < 3\varepsilon ] > 1-\delta$.
:::


:::{prf:proof}
:class: dropdown
:enumerated: false
Let $G = \{ i : \| \vec{x}^* - \vec{x}_i \| \leq \varepsilon\}$.
Since $\mathbb{E}[|G|] > 2 m/3$, it follows from a standard [Chernoff Bound](https://en.wikipedia.org/wiki/Chernoff_bound#Multiplicative_form_(relative_error)) that $\Pr(|G| \leq t/2)\leq \delta$ for $m = O(\log(1/\delta))$.
Therefore, it suffices to show that if $|G|> m/2$, then $\| \vec{x}^*-\widehat{\vec{x}} \| < 3\varepsilon$.

Suppose $|G| > m/2$.
By the triangle inequality, for each $i,j\in[m]$,
$$
    d(i,j) = \|\vec{x}_i - \vec{x}_j\|
    \leq \| \vec{x}^* - \vec{x}_i \| + \| \vec{x}^* - \vec{x}_j \|.
$$
Therefore, for each $i\in G$
$$
    \big|\{j\in [m] : d(i,j) \leq 2\varepsilon\}\big|
    > t/2,
$$
and hence, 
\[
c(i) = \operatorname{median}\big(d(i,1),\ldots, d(i,t)\big) \leq 2\varepsilon.
\]

It follows that $c(i^*) = \min_i c(i) \leq 2\varepsilon$, and since $c(i^*) = \operatorname{median}(d(i^*,1),\ldots, d(i^*,t))$, we have
$$
    \big|\{j\in [m] : d(i^*,j) \leq 2\varepsilon\}\big|
    \geq t/2.
$$
But also $|G|>t/2$. 
So, by the pigeonhole principle, there is at lest one index $j^*\in G$ for which 
$$
    d(i^*,j^*) \leq 2\varepsilon.
$$
Therefore, by the triangle inequality
$$
    \| \vec{x}^* - \vec{x}_{i^*} \|
    \leq \| \vec{x}^* - \vec{x}_{j^*} \| + \| \vec{x}_{i^*} - \vec{x}_{j^*} \|
    \leq \varepsilon + 2\varepsilon = 3\varepsilon,
$$
as desired.
:::



