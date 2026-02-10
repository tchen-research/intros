---
title: Quantum Typicality
description: What is quantum typicality?
bibliography: typicality.bib
link-citations: true
---


## Introduction


Randomized algorithms for approximating the trace of a $n\times n$ matrix $\vec{A}$ using only matrix-vector products with $\vec{A}$ are increasingly used as a computational primitive in a number of areas in the computational sciences.


Many commonly used stochastic trace estimation algorithms make critical use of the fact that $\bm{\varphi}^\T \vec{A} \bm{\varphi}$ is an unbiased estimator for $\tr(\vec{A})$, at least provided that $\EE[\bm{\varphi} \bm{\varphi}^\T] = \vec{I}$.
Indeed,
$$
    \EE[ \bm{\varphi}^\T\vec{A}\bm{\varphi} ] 
    = \EE[ \tr(\bm{\varphi}^\T\vec{A}\bm{\varphi}) ] 
    = \EE[ \tr(\vec{A}\bm{\varphi}\bm{\varphi}^\T) ] 
    =  \tr(\vec{A}\EE[\bm{\varphi}\bm{\varphi}^\T])
    = \tr(\vec{A}).
$$
Common choice of distribution for $\bm{\varphi}$ include sampling the entries independently from $\operatorname{Unif}(\{-1,+1\})$ or $\mathcal{N}(0,1)$ or sampling $\bm{\varphi}$ from $\sqrt{n}\operatorname{Unif}(\mathbb{S}^{n-1})$, the uniform distribution on the hypersphere of radius $\sqrt{n}$.
For all of these distributions, tail bounds of the form 
$$
    \PP\Big[\big|\bm{\varphi}^\T \vec{A} \bm{\varphi} - \tr(\vec{A})\big| > \epsilon\Big] 
    \leq \delta = \delta(\vec{A}, \epsilon)
$$

have been studied [@reimann_07;@popescu_short_winter_06;@avron_toledo_11;@roostakhorasani_ascher_14;@meyer_musco_musco_woodruff_21;@cortinovis_kressner_21;@chen_trogdon_ubaru_22].
Current state of the art bounds offer refined tail bounds for $\delta$ depending on $\vec{A}$ through $\|\vec{A}\|_\F^2$ and $\|\vec{A}\|_2$ with the precise bounds depending on the choice of distribution for $\vec{v}$.

In this note, we provide a brief historical overview of these algorithms and their corresponding analyses. 
While our exposition is far from comprehensive, we believe it provides some insight into the history of typicality-based algorithms.
We hope that some of the connections highlighted in this letter motivate further the interaction between physics, numerical analysis, and theoretical computer science.

## Quantum mechanics and linear algebra


In the now standard mathematical formalism of quantum mechanics, pure states of a quantum system are represented as (unit length) elements of a Hilbert space $\mathcal{H}$, and physical observables (such as energy, translational momentum, orbital angular momentum, and spin) are self-adjoint operators on $\mathcal{H}$.
If $\mathcal{H}$ has dimension $n<\infty$ (which we shall assume throughout this paper), the spectral theorem ensures an observable $\vec{A}$ can be decomposed as
$$
 \vec{A} = \sum_{i=1}^{n} \lambda_i \bm{\psi}_i \bm{\psi}_i^\T,
$$
where $\{\lambda_i\}$ are the eigenvalues of $\vec{A}$ and $\{ \bm{\psi}_i \}$ are the (orthonormal) eigenstates.

When a measurement is made of $\vec{A}$ while the system is in a state $\bm{\psi}_i$, the value $m(\vec{A};\bm{\psi}_i)$ of the measurement observed will be $\lambda_i$.
More generally, if the measurement is made while the system is in an arbitrary state $\bm{\varphi}$, the value of the measurement observed may be any of $\lambda_1, \ldots, \lambda_n$.
In particular, the values observed will be $\lambda$ with probability proportional to the norm of the projection of $\bm{\varphi}$ onto the the eigenstates associated with $\lambda$.
That is, 
$$
    \PP_{\textup{qm}}\big[ m(\vec{A};\bm{\varphi}) = \lambda \big] = \sum_{i:\lambda_i = \lambda} |\bm{\psi}_i^\T\bm{\varphi}|^2.
    $$
Here the subscript ``qm'' indicates that the probability is with respect to the inherent randomness in the formalism of quantum mechanics.

That the results of a measurement are probabilistic is distinctly quantum; in the classical setting, repeated measurements of a system in a given state result in exactly the same measurement value.
However, repeated measurements of an observable in a given state will eventually reveal an average value
$$
\langle \vec{A} \rangle_{\bm{\varphi}}
    := \EE_{\textup{qm}}[m(\vec{A};\bm{\varphi})]
    = \sum_{i=1}^{n} \lambda_i |\bm{\psi}_i^\T\bm{\varphi}|^2
    = \bm{\varphi}^\T \vec{A} \bm{\varphi},
$$
referred to as \emph{quantum expectation value} (QEV) of $\vec{A}$ in state $\bm{\varphi}$.

In fact, as noted by von Neumann [@vonneumann_29]
, the probability distribution \cref{eqn:Aphi_dist} can be recovered by knowledge of the moments $\langle \vec{A} \rangle_{\bm{\varphi}}, \langle \vec{A}^2 \rangle_{\bm{\varphi}}, \langle \vec{A}^3 \rangle_{\bm{\varphi}}, \ldots$

> The simplest description of a state by means of a wave function $\varphi$ is obtained in this way: the expectation value of the quantity $A$ in the state $\varphi$ is equal to $(A\varphi, \varphi)$. The specification of all expectation values provides, as it includes the expectation values of all powers (i.e., the so-called higher moments of a probability distribution), knowledge of the entire probability distribution of every quantity—and thus a complete statistical characterization of the system

More generally, a quantum system may be be in a mixture of pure states $\{\bm{\varphi}_i\}$ each occurring with probability $\{p_i\}$.
Such an ensemble is represented mathematically by the density matrix
$$
    \vec{\rho} = \sum_{i=1}^{n} p_i \bm{\varphi}_i \bm{\varphi}_i^\T.
$$
When the observable $\vec{A}$ is measured in such a state, the corresponding QEV is given by
$$
    \langle \vec{A} \rangle_{\vec{\rho}}
    %:= \EE[m(\vec{A},\vec{\rho})]
    = \sum_{i} \bm{\varphi}_i^\T \vec{A} \bm{\varphi}_i p_i
    = \tr(\vec{\rho} \vec{A}).
$$
Note that in the case of a single pure state, so that $\vec{\rho} = \bm{\varphi}\bm{\varphi}^\T$, we recover the original formula for the QEV by the identity $\tr(\bm{\varphi}\bm{\varphi}^\T\vec{A}) = \bm{\varphi}^\T \vec{A} \bm{\varphi}$.
