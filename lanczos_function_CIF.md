---
title: Error bounds for Lanczos-based matrix function approximation
description: We provide some new error bounds for the Lanczos methods for computing matrix functions.
bibliography: lanczos_function_CIF.bib
link-citations: true
math:
  '\lan': '\mathsf{lan}'
  '\err': '\mathsf{err}'
  '\Res': '\mathsf{res}'
---


This is a companion piece to @chen_greenbaum_musco_musco_22. 
This approach has prompted a number of follow-up works that generalize it to non-symmetric problems, and block and rational algorithms @xu_chen_24, @simunec_24, @adler_hu_pan_xue_25.

Code to replicate the figures in the corresponding paper can be found on [Github](https://github.com/tchen-research/lanczos_function_CIF).

## Introduction

The Lanczos method for matrix function approximation ([Lanczos-FA](./lanczos-fa.html)) can be used to approximate $f(\vec{A})\vec{b}$, and in the case that $f(x) = 1/x$ and $\vec{A}$ is positive definite, this approximation is optimal over Krylov subspace.
This case is very well studied and a range of error bounds and estimates exist.
However, for other functions, the standard bounds are often too pessimistic as they do not take into account fine grained information about the spectrum of $\vec{A}$ such as outlaying or clustered eigenvalues.
This makes it difficult to know when Lanczos-FA has reached a suitable accuracy for a given problem. 


In this paper we show how to reduce the error of approximating $f(\vec{A})\vec{b}$ with Lanczos-FA to the error of solving a certain linear system with the Lanczos-FA.
This allows us to leverage the range of existing bounds for the convergence of Lanczos-FA on linear systems to easily obtain a priori and a posteriori bounds for other matrix functions including piecewise analytic functions such as the sign function. 
Our a posteriori bounds are highly accurate and can be used as practical stopping criteria. 

## The basic idea

The $k$-th Lanczos-FA approximation to $f(\vec{A}) \vec{b}$ is defined as
$$
    \lan_k(f)
    := \vec{Q} f(\vec{T}) \vec{Q}^{\T} \vec{b},
$$
where $\vec{Q}$ and $\vec{T}$ are produced by the Lanczos method run for $k$ steps on $(\vec{A},\vec{b})$.

Suppose $\vec{A}$ is a Hermitian matrix and If $f$ is analytic in a neighborhood of the eigenvalues of $\vec{A}$ and $\Gamma$ is a contour in this neighborhood containing these eigenvalues,
$$
    f(\vec{A})\vec{b} = - \frac{1}{2 \pi i} \oint_{\Gamma} f(z) (\vec{A} - z \vec{I} )^{-1} \vec{b} \, \d{z}.
$$
If the eigenvalues of $\vec{T}$ are also contained in $\Gamma$, we similarly have
$$
    \vec{Q} f(\vec{T}) \vec{Q}^\T \vec{b}
    = - \frac{1}{2 \pi i} \oint_{\Gamma} f(z) \vec{Q} (\vec{T} - z \vec{I} )^{-1} \vec{Q}^\T \vec{b} \, \d{z} .
$$
Combining, these, we see that the Lanczos-FA error can be written as
$$
f(\vec{A}) \vec{b} - \vec{Q} f( \vec{T} ) \vec{Q}^{\T} \vec{b}
    = - \frac{1}{2 \pi i} \oint_{\Gamma} f(z) \, \err_k(z) \, \d{z}.
$$

For $z \in \mathbb{C}$, define the $k$-th Lanczos-FA error and residual for the linear system $(\vec{A}-z\vec{I}) \vec{x} = \vec{b}$ as,
$$\begin{align*}
    \err_k(z,\vec{A},\vec{b}) &:= (\vec{A} - z \vec{I})^{-1} \vec{b} - \vec{Q}(\vec{T}-z\vec{I})^{-1}\vec{Q}^\T \vec{b}%\lan_k ( h_z )
    ,\\
    \Res_k(z,\vec{A},\vec{b}) &:= \vec{b} - (\vec{A} - z \vec{I}) \vec{Q}(\vec{T}-z\vec{I})^{-1}\vec{Q}^\T \vec{b}.%\,\lan_k ( h_z ).
\end{align*}$$

It is a well-known fact that
$$
    \Res_k(z) 
    = \left( \frac{(-1)^{k}}{\det(\vec{T} -z \vec{I}) }\prod_{j=1}^{k} \beta_j \right) \| \vec{b} \|_2\: \vec{q}_{k+1}.
$$
Consequently, for all $z , w \in \mathbb{C}$, where $\vec{A} - z \vec{I}$ and $\vec{A} - w \vec{I}$ are both invertible,
$$\begin{align*}
    \err_k(z)
    &= \det(h_{w,z}(\vec{T}))  h_{w,z}(\vec{A}) 
    \,\err_k(w)
    \\
    \Res_k(z)
    &= \det(h_{w,z}(\vec{T}))
    \,\Res_k(w),
\end{align*}$$
where 
$$
h_{w,z}(x) = \frac{x-w}{x-z}.
$$

Thus, if $\Gamma$ is a simple closed curve or union of simple closed curves inside this neighborhood and enclosing the eigenvalues of $\vec{A}$ and $\vec{T}$ and $w$ a point not in $\Lambda(\vec{T})\cup\Lambda(\vec{A})$,
$$
    f(\vec{A}) \vec{b} - \lan_k(f)
    = \left( - \frac{1}{2\pi i} \oint_{\Gamma} f(z) \det(h_{w,z}(\vec{T})) h_{w,z}(\vec{A}) \d{z} \right) \, \err_k(w).
$$

Applying the triangle inequality and the submultiplicitivity of matrix-norms, we can then obtain a bound  
$$
    \| f(\vec{A})\vec{b} - \lan_k(f) \|
    \leq \underbrace{\vphantom{ \bigg| }\left( \frac{1}{2\pi}\oint_{\Gamma} |f(z)| \left(\prod_{i=1}^{k} \| h_{w,z}\|_{S_i}\right) \|h_{w,z}\|_{S_0} | \d{z} | \right)}_{\text{integral term}}  \hspace{-.5 em}\underbrace{\vphantom{ \Bigg| } \| \err_k(w) \| , \hspace{-.4em} }_{\text{linear system error}} \hspace{-.5em}
$$
where $S_i$ are some suitably chosen sets and $\|g\|_S:= \max_{x\in S}|g(x)|$.

Note that the integral term and linear system error term in the theorem are entirely decoupled!
Thus, once the integral term is computed, bounding the error of Lanczos-FA for $f(\vec{A})\vec{b}$ is reduced to bounding $\| \err_k(w) \|$, and if the integral term can be bounded independently of $k$, This implies that, up to a constant factor, the Lanczos-FA approximation to $f(\vec{A})\vec{b}$ converges at least as fast as $\| \err_k(w) \|$.

## What else is in the paper?

Much of the paper is focused on practical aspects regarding the use of such a bound.
In particular:

- we provide a detailed discussion and analysis of the validity of our bounds when the Lanczos iterate was computed using finite precision arithmetic.
- we provide several analytic examples where the integral term is computed directly
- we use numerical experiments to demonstrate the effectiveness of our bounds on a variety of functions such as the square root, inverse square root, and sign functions, and
- we derive similar bounds for quadratic forms.

It's worth noting that similar ideas have been previously used to get error bounds for Stieltjes functions [@ilic_turner_simpson_09;@frommer_guttel_schweitzer_14;@frommer_schweitzer_15]. 
As discussed in Section 2.2 of the paper, our bounds differ in a number of key ways.
One key difference is that we reduce error bounds for Lanczos-FA a general function to error bounds for Lanczos-FA on a fixed linear system.
This allows intuition about the [convergence of Lanczos-FA on linear systems](./lanczos-fa_linear_systems.html) to be transferred to other functions.
In addition, our analysis allows bounds for piecewise continuous functions like the sign function.
