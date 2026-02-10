---
title: Near-optimal convergence of the full orthogonalization method
description: We establish a near-optimality guarantee for the full orthogonalization method (FOM), showing that the overall convergence of FOM is nearly as good as GMRES. 

bibliography: FOM_convergence.bib
math: 

---


This is a companion piece to @chen_meurant_24.

Code to replicate the figures in the corresponding paper can be found on [Github](https://github.com/tchen-research/FOM_convergence).


## Introduction   


The *full orthogonalization method* (FOM) @saad_81 and the *generalized minimal residual method* (GMRES) @saad_schultz_86 are two Krylov subspace methods used for solving a real or complex non-symmetric  linear system of equations.
Note that if $\vec{A}$ is symmetric GMRES is mathematically equivalent to MINRES @paige_saunders_75, and if $\vec{A}$ is symmetric positive definite FOM is mathematically equivalent to conjugate gradient @hestenes_stiefel_52. While this is relevant for efficient implementation, it does not impact the exact arithmetic theory in this note.
\begin{equation}
    \label{eqn:Axb}
\vec{A}\vec{x} = \vec{b}.
\end{equation}
Assuming an initial guess $\vec{x}_0 = \vec{0}$, both FOM and GMRES produce iterates from the Krylov subspace
\begin{equation}
\label{eqn:krylov_subspace}
    \mathcal{K}_k(\vec{A},\vec{b}) 
    := \operatorname{span}\{\vec{b}, \vec{A}\vec{b}, \ldots, \vec{A}^{k-1}\vec{b}\}.
\end{equation}
but according to slightly different formulas.


It is well-known that the GMRES iterates satisfy a residual optimality guarantee:
\begin{equation}\label{eqn:GMRES_opt}
    \vec{x}_k^{\mathrm{G}} 
    = \operatornamewithlimits{argmin}_{\vec{x}\in\mathcal{K}_k(\vec{A},\vec{b})} \| \vec{b} - \vec{A} \vec{x} \|_{2}.
\end{equation}
Hence, the GMRES residual norms are non-increasing and are optimal among Krylov subspace methods.
This optimality guarantee leads to well-known results on the rate of convergence of the residual norms in terms of quantities such as the condition number of $\vec{A}$ @saad_03.
On the other hand, the FOM residual norms often appear oscillatory, with large jumps.
In fact, it is easy to construct examples for which $\|\vec{r}_k^{\mathrm{F}}\|_2/\|\vec{r}_0^{\mathrm{F}}\|_2$ can be arbitrarily large @meurant_tebbens_20!

The main result of our paper is to prove that, in an *overall sense*, FOM behaves nearly as well as GMRES:

:::{prf:theorem}
For every $k\geq 1$,
\begin{equation*}
    \min_{0\leq j\leq k} \|\vec{r}_{j}^{\mathrm{F}}\|_2
    \leq \sqrt{k+1} \cdot \|\vec{r}_k^{\mathrm{G}}\|_2.
\end{equation*}
:::

The paper is actually very short and easy to read, so I'm not including that much more here.