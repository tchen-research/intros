---
title: Numerical computation of the equilibrium-reduced density matrix for strongly coupled open quantum systems
description: We describe a numerical algorithm to approximate partial traces. This can be used to approximate reduced density matrices of thermal states.
bibliography: mean_force.bib
math: 
  '\sysa': '\mathrm{a}'
  '\sysb': '\mathrm{b}'
  '\syst': '\mathrm{t}'
  '\da': 'd_{\sysa}'
  '\db': 'd_{\sysb}'
  '\dt': 'd_{\syst}'
  '\Ha': '\mathcal{H}_{\sysa}'
  '\Hb': '\mathcal{H}_{\sysb}'
  '\Ht': '\mathcal{H}_{\syst}'
  '\trb': '\tr_{\sysb}'
  '\trbest': '\widehat{\mathsf{tr}}_{\sysb}'
---


This is a companion piece to @chen_cheng_22. See also @chen_chen_li_nzeuton_pan_wang_24 for some follow-up work.

Code to replicate the figures in the corresponding paper can be found on [Github](https://github.com/tchen-research/mean_force).


## Introduction   

A quantum system is described by a Hamiltonian matrix $\vec{H}$.
When the system is at thermal equilibrium at temperature $1/\beta$, then the state of the system is described by the matrix
$$
\bm{\rho} = \exp(-\beta \vec{H}) / Z(\beta)
,\qquad 
Z(\beta) = \tr(\exp(-\beta \vec{H})).
$$
The function $Z(\beta)$ is called the partition function, and gives us access to all kinds of thermodynamic quantities such as the energy, specific heat, magnetization, etc.

Often we are interested in the state of some subsystem of the total system. 
This is represented by 
$$
\bm{\rho}^*(\beta) = \trb(-\beta \vec{H}) / Z(\beta),
$$
where $\trb(\cdot)$ is the *partial trace* described below.


### Exponentially big linear algebra

Computing the trace (and partial trace) of a given matrix is trivial. For the trace you just sum the diagonal entries, and for the partial trace you do something similarly easily.
However, we are starting with $\vec{H}$ not $\exp(-\beta \vec{H})$, and this presents some issues.

For a system with $N$ sites, the size of $\vec{H}$ is $2^N$. 
This means that even for reasonably small $N$, storing a $N\times N$ matrix is impossible.
For instance, if $N=20$, then $2^N$ is about one million, and storing a $N\times N$ matrix of double precision floating point numbers (64 bits per number) would require over *8 terrabtyes* of memory!
Increasing $N$ by one quadruples the memory costs, so even if we have 8 terrabytes of memory, we'll quickly run out for $N$ even just a bit larger.

Fortunately, $\vec{H}$ typically has a very sparse representation, often with just $O(N)$ nonzero entries.
While this means we can store $\vec{H}$ itself, $\exp(-\beta \vec{H})$ is typically not sparse and therefore it is intractable to store.
Moreover, even if we would store $\exp(-\beta \vec{H})$, computing it would be a big task!


### Stochastic trace estimation

There are lots of algorithms for estimating $\vec{v}^\T f(\vec{H})\vec{v}$ using only matrix-vector products with $\vec{H}$. 
For instance, the Lanczos method for matrix function approximation ([Lanczos-FA](./lancos-fa.html)) is a common choice.

There are lots of nice overviews on stochastic trace estimation, and this [blog post](https://www.ethanepperly.com/index.php/2023/01/26/stochastic-trace-estimation/) by Ethan Epperly is especially good.
For our purposes, we only need to know that if $\vec{v}$ is a vector for which $\EE[\vec{v} \vec{v}^\T] = \vec{I}$ (i.e. all the entries are uncorrelated), then
$$
\EE\big[\vec{v}^\T \vec{A} \vec{v}\big] = \tr(\vec{A}).
$$
A simple way to generate such a $\vec{v}$ is by taking each entry to be an independent standard normal random variable.
We will call estimators of the form $\vec{v}^\T \vec{A}\vec{v}$ a *quadratic trace estimator*.


## Partial trace estimation

Suppose we can partition a matrix $\vec{A}$ as
$$
    \vec{A}
    = 
    \begin{bmatrix}
    \vec{A}_{1,1} & \vec{A}_{1,2} & \cdots & \vec{A}_{1,\da} \\
    \vec{A}_{2,1} & \vec{A}_{2,2} & \cdots & \vec{A}_{2,\da} \\
    \vdots & \vdots & \ddots & \vdots \\
    \vec{A}_{\da,1} & \vec{A}_{\da,2} & \cdots & \vec{A}_{\da,\da}
    \end{bmatrix}.
$$
Then the partial trace of $\vec{A}$ is defined as
$$
    \trb(\vec{A}) :=
    \begin{bmatrix}
    \tr(\vec{A}_{1,1}) & \tr(\vec{A}_{1,2}) & \cdots & \tr(\vec{A}_{1,\da}) \\
    \tr(\vec{A}_{2,1}) & \tr(\vec{A}_{2,2}) & \cdots & \tr(\vec{A}_{2,\da}) \\
    \vdots & \vdots & \ddots & \vdots \\
    \tr(\vec{A}_{\da,1}) & \tr(\vec{A}_{\da,2}) & \cdots & \tr(\vec{A}_{\da,\da})
    \end{bmatrix}.
$$

Note that each entry of the partial trace is obtained by taking the trace of a certain matrix.
This suggests we can apply use quadratic trace estimators to approximate each of the entries of the partial trace!
Indeed, 
$$
 (\vec{I}_{\da} \! \otimes \vec{v})^\T
    \vec{A}
    (\vec{I}_{\da} \! \otimes \vec{v})
    =     
    \begin{bmatrix}
    \vec{v}^\T \vec{A}_{1,1} \vec{v} & \vec{v}^\T \vec{A}_{1,2} \vec{v} & \cdots & \vec{v}^\T \vec{A}_{1,\da} \vec{v} \\
    \vec{v}^\T \vec{A}_{2,1} \vec{v} & \vec{v}^\T \vec{A}_{2,2} \vec{v} & \cdots & \vec{v}^\T \vec{A}_{2,\da} \vec{v} \\
    \vdots & \vdots & \ddots & \vdots \\
    \vec{v}^\T \vec{A}_{\da,1} \vec{v} & \vec{v}^\T \vec{A}_{\da,2} \vec{v} & \cdots & \vec{v}^\T \vec{A}_{\da,\da} \vec{v}
    \end{bmatrix}
    $$
provides an unbiased estimator for $\trb(\vec{A})$.

Given independent and identically distributed copies $\vec{v}_1, \ldots, \vec{v}_m$ of $\vec{v}$, we arrive at an estimator
$$
    \trbest^{m}(\vec{A})
    := \frac{1}{m} 
    \sum_{i=1}^{m}
    (\vec{I}_{\da} \! \otimes \vec{v}_i)^\T \vec{A} (\vec{I}_{\da} \! \otimes \vec{v}_i).
$$



