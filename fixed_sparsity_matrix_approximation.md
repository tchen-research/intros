---
title: Fixed-sparsity matrix approximation from matrix-vector products
description: We provide some new error bounds for the Lanczos methods for computing matrix functions.
bibliography: fixed_sparsity_matrix_approximation.bib
link-citations: true
---


This is a companion piece to @amsel_chen_halikias_dumankeles_musco_musco_26.

Code to replicate the figures in the corresponding paper can be found on [Github](https://github.com/tchen-research/fixed_sparsity_matrix_approximation).



## Introduction

In this work, we study the task of approximating a $d\times d$ matrix $\vec{A}$ with a matrix $\widetilde{\vec{A}}$ of a *pre-specified* sparsity pattern (e.g. diagonal, tridiagonal, block diagonal, etc.).
There are many existing linear algebra algorithms for matrices of a given sparsity pattern, so this task yields a proxy for $\vec{A}$ which is compatible with such algorithms. 
For example, a banded approximation to $\vec{A} = \vec{M}^{-1}$ could be used as a preconditioned to systems involving $\vec{M}$.


The catch is that we assume $\vec{A}$ can only be accessed using matrix-vector products (called *matvecs*).
That is, we restrict ourselves to algorithms which choose vectors $\vec{x}_1, \ldots, \vec{x}_m$ and subsequently receive responses $\vec{A}\vec{x}_1, \ldots, \vec{A}\vec{x}_m$ from some black-box, but cannot look at $\vec{A}$ in any other way.
We will measure the cost of algorithms by $m$, the total number of matrix-vector queries done, which in many settings may be a realistic proxy for the real-world costs.
This is often called the *matrix-vector query model*.

It's worth noting that we can recover $\vec{A}$ in this access model by simply reading off the columns of $\vec{A}$ one-by-one.
This means we can effectively run any classical algorithm for the cost of $d$ matvecs.
If $d$ is very big this is not practical, but it is worth keeping in mind as a trivial upper bound for the cost of any algorithm in this model.
In general, our goal will be to design algorithms which use $m\ll d$ matvecs.


One of the most powerful features of the matvec query model is that it is often possible to prove *lower-bounds* for the hardness of a given task.
That is, to show that it is impossible (for *any* algorithm) to solve the task unless the algorithm uses at least a certain number of queries.
This is exciting, because it offers a way of understanding the "difficulty" of linear algebra tasks!
In contrast, proving lower-bounds on the number of floating point operations (which has classically been the measure of cost for linear-algebra algorithms), is seemingly far beyond current techniques. 


One thing I think is really cool about this paper, is that almost everything is self-contained!
In particular, besides a few well-known results on properties of Gaussian matrices and a qualitative central limit theorem called the Berry-Esseen theorem, everything needed is proved within the paper using relatively elementary techniques which would be understandable to anyone with an undergraduate-level knowledge of linear-algebra and probability.
In fact, the general technique is so simple that we even prove a lower-bound for diagonal estimation on this page!
Overall, I hope that this paper can serve as a case-study on the power and beauty of the matrix-vector query model.


## Two theoretical problems

For convenience, we will measure the error in the Frobenius norm:
$$
\|\vec{A} - \widetilde{\vec{A}} \|_{\F}
= \sqrt{\sum_{i,j} \big([\vec{A}]_{\,i,j} - [\widetilde{\vec{A}}]_{\,i,j} \big)^2}.
$$
Let's also encode the desired sparsity pattern in a binary matrix $\vec{S}$.
The nonzero entries of $\vec{S}$ will tell us where $\widetilde{\vec{A}}$ is allowed to have a nonzero entries.
It is pretty clear then that the best sparse approximation to $\vec{A}$ is simply $\vec{S}\circ \vec{A}$, where "$\circ$" is the Hadamard (entrywise) product. 
That is, the best approximation is obtained by zeroing out the entries of $\vec{A}$ which are required to be zero.


We will study the following problems:

:::{prf:definition} Best approximation by a matrix of fixed sparsity
Given a matrix $\vec{A} \in\mathbb{R}^{n\times d}$ and a binary matrix $\vec{S} \in \{0,1\}^{n\times d}$, find a matrix $\widetilde{\vec{A}}$ so that $\widetilde{\vec{A}} = \vec{S}\circ\widetilde{\vec{A}}$ and
$$
    \| \vec{A} - \widetilde{\vec{A}} \|_\F
    \leq (1+\varepsilon) \, \| \vec{A} - \vec{S}\circ \vec{A} \|_\F.
$$
:::

:::{prf:definition} Best approximation to on-sparsity-pattern entries
Given a matrix $\vec{A} \in\mathbb{R}^{n\times d}$ and a binary matrix $\vec{S} \in \{0,1\}^{n\times d}$, find a matrix $\widetilde{\vec{A}}$ so that $\widetilde{\vec{A}} = \vec{S}\circ\widetilde{\vec{A}}$ and
$$
    \| \vec{S}\circ \vec{A} - \widetilde{\vec{A}} \|_\F^2
    \leq \varepsilon \,\| \vec{A} - \vec{S}\circ \vec{A} \|_\F^2.
$$
:::

Note that since $\widetilde{\vec{A}}$ has the same sparsity as $\vec{S}$; i.e. $\widetilde{\vec{A}} = \vec{S}\circ\widetilde{\vec{A}}$,
$$
\| \vec{A} - \widetilde{\vec{A}} \|_\F^2
= \| \vec{A} - \vec{S}\circ \vec{A} \|_\F^2 + \|\vec{S}\circ\vec{A} - \widetilde{\vec{A}}\|_\F^2.
$$
This allows us to relate the two problems, and for $\varepsilon \in (0,1)$, the two are equivalent up to (small) constants.


## A simple algorithm

When the sparsity pattern $\vec S$ has at most $s$ non-zero entries per row, this algorithm uses $m = O(s/\varepsilon)$ non-adaptive matrix-vector product queries.
Specifically, the algorithm computes $\vec{Z} = \vec{A}\vec{G}$, where $\vec{G}$ is a $d\times m$ matrix with independent standard normal entries, and then outputs the matrix 
$$
    \widetilde{\vec{A}} = \operatornamewithlimits{argmin}_{\vec{X} = \vec{S}\circ\vec{X}} \| \vec{Z} - \vec{X} \vec{G}\|_\F.
$$
It's not hard to see how to implement this algorithm. 
First, observe that each row of the approximation can be computed independently. 
Then, note that the solution to each row is obtained by solving a standard least squares problem involving the relevant rows of $\vec{G}$.

Using standard tools from high dimensional probability, which characterize the expected Frobenius norm of the pseudoinverse of a Gaussian matrix, it is relatively simple to establish the following convergence guarantee:

:::{prf:theorem}
:label: thm:upper-bound
Consider any $\vec{A} \in\mathbb{R}^{n\times d}$ and any $\vec{S}\in\{0,1\}^{n\times d}$ with at most $s$ nonzero entries per row. 
Then, for any $m\geq s+2$, using $m$ randomized matrix-vector queries, the matrix $\widetilde{\vec{A}}$ is equal to $\vec{S}\circ\vec{A}$ in expectation and satisfies
$$
\EE\Bigl[\|\vec{S}\circ\vec{A} - \widetilde{\vec{A}}\|_\F^2 \Bigr]
\leq {\frac{s}{m-s-1}} \| \vec{A} - \vec{S}\circ \vec{A} \|_\F^2.
$$
The above inequality is equality if each row of $\vec{S}$ has exactly $s$ non-zero entries.
:::

If we set $m = O(s/\varepsilon)$, then we solve the two problems above, at least in expectation. 
One can then use Markov's inequality to get a probability bound.
This matches the gurantees for a number of existing algorithms for special cases of Problem 1 and 2 [@curtis_powell_reid_74,], [@bekas_kokiopoulou_saad_07], [@halikias_townsend_23], [@dharangutte_musco_23], etc., and we comment on the connections in more detail in the paper.


## A lower bound

For any fixed $\vec{A}$, there is an algorithm which solves Problem 1 exactly using zero queries; simply hard-code the algorithm to output $\vec{S}\circ\vec{A}$.
This is not a practically useful algorithm, because it only works for the one particular $\vec{A}$ and we would be unable to practically ``find'' such an algorithm for a given matrix.
However, it means it's impossible to prove a statement like: "For any matrix $\vec{A}$, no algorithm using (some amount of queries depending on $\varepsilon$) can solve Problem 1".
Instead, it is standard to consider a *distribution* of matrices $\vec{A}$ (or equivalently, a random matrix). 
While an algorithm might be hard-coded for the particular distribution we consider, there is hope that the distribution is always still very random, even conditioned on a small number of arbitrary matvec queries.


In particular, in the paper we prove a precise version of the following:

:::{prf:theorem}
:label: thm-lower
There are constants $c,C>0$ so that, for any $\varepsilon\in(0,c)$ and $s\geq 1$, there is a distribution on matrices $\vec{A}$ such that, for any sparsity pattern $\vec{S}$ in a broad class of sparsity patterns, and for any algorithm that uses $m\leq Cs/\varepsilon$ matrix-vector queries to $\vec{A}$ to output a matrix $\widetilde{\vec{A}}$, 
$$
    \PP\Bigl[ \| \vec{A} - \widetilde{\vec{A}} \|_\F
    \leq (1+\varepsilon) \,\| \vec{A} - \vec{S}\circ \vec{A} \|_\F \Bigr] \leq \frac{1}{25}.
$$
:::

This shows that there are problems for which no algorithm can use more than a constant (e.g. $100$) times fewer matvecs than the algorithm analyzed in the previous section.

At a high-level, the lower-bound technique is pretty simple. 
The basic idea is that if $\vec{A} = \vec{G}^\T\vec{G}$, where $\vec{G}$ is $d\times d$ random Gaussian matrix, then after a small number of (possibly adaptive) matvec queries to $\vec{A}$, there is still too much randomness in $\vec{A}$ (conditioned on the information leared by the queri– in the compressed sensing setting, we have a necessaryes) to learn what we need [@braverman_hazan_simchowitz_woodworth_20].

## A proof for diagonal estimation

To give a flavor of how the above lower-bound above is proved, we will now prove a simpler result for diagonal estimation; i.e. for when $\vec{S} = \vec{I}$.
We will be a bit cavalier 🤠 , providing a TCS style proof (you shouldn't be just trusting things on the internet anyway).

:::{prf:theorem}
:label: thm:lower-diag
There is a constant $C>0$ such that, for any $\varepsilon>0$ there are matrices $\vec{A}$ such that, for any algorithm that uses $m\leq C/\varepsilon$ matrix-vector queries to $\vec{A}$ to output a diagonal matrix $\widetilde{\vec{A}}$, 
$$
    \PP\Bigl[ \| \vec{I}\circ \vec{A} - \widetilde{\vec{A}} \|_\F^2
    \leq \varepsilon \,\| \vec{A} - \vec{I}\circ \vec{A} \|_\F^2 \Bigr] \leq \frac{1}{25}.
$$
:::

The key technical tool will be the following Lemma:

:::{prf:lemma}
Suppose $\vec{A}\sim\operatorname{Gaussian}(d,d)$ and let $\vec{x}_1, \ldots, \vec{x}_m$ and $\vec{y}_1 = \vec{A}\vec{x}_1, \ldots, \vec{y}_m = \vec{A}\vec{x}_m$ be such that, for each $j=0,1,\ldots, m$, $\vec{x}_j$ was chosen based only on the query vectors $\vec{x}_1, \ldots, \vec{x}_{j-1}$ and the outputs $\vec{y}_1, \ldots, \vec{y}_{j-1}$ and (w.l.o.g.) unit-length and orthogonal to $\vec{x}_1, \ldots, \vec{x}_{j-1}$.

Let $\vec{X} = [\vec{x}_1, \ldots, \vec{x}_m]$ and $\vec{Y} = [\vec{y}_1, \ldots, \vec{y}_m]$. 
Then, $\vec{A}$ can be factored as
$$
\vec{A} = \begin{bmatrix}\vec{Y} & \vec{G} \end{bmatrix} \begin{bmatrix}\vec{X}^\T \\ \vec{Z}^\T \end{bmatrix},
$$
where $\vec{Z}^\T\vec{X} = \vec{0}$, $\vec{Z}^\T\vec{Z} = \vec{I}$, and $\vec{G}\sim\operatorname{Gaussian}(d,d-m)$, independently of $\vec{X}$ and $\vec{Y}$.
:::

:::{prf:proof}
:class: dropdown
:enumerated: false
Clearly $\vec{A}$ must be decomposed as above for \emph{some} matrix $\vec{G}$ (because this is just encoding what $\vec{A}$ does to the vectors in $\vec{X}$).
It therefore suffices to argue that $\vec{G}$ is a Gauassian matrix independent of $\vec{X}$ and $\vec{Y}$.

We proceed by induction, noting that the base case is trivial. 

Suppose the lemma holds for $m$ queries.
Choose $\vec{x}_{m+1}$ as a unit vector orthogonal to $\vec{x}_1, \ldots, \vec{x}_m$ (i.e. so that $\vec{X}^\T\vec{x}_{m+1} = \vec{0}$.
Note that $\vec{x}_{m+1} = \vec{Z} \vec{c}$, for some length $d-m$ unit vector $\vec{c}$; indeed, if $\vec{x}_{m+1}$ is orthogonal to $\vec{X}$ then it is in the span of $\vec{Z}$, and it is unit-length if and only if $\vec{c}$ is unit-length.


Deterministically complete $\vec{c}$ or an orthogonal matrix\footnote{For instance, append the identity, delete the first linearly dependent column, and apply Gram--Schmidt.}
$$
\vec{C} = \begin{bmatrix}
    \vec{c} & \hat{\vec{C}}
\end{bmatrix}.
$$
By the invariance of Gaussian vectors under unitary transforms, 
$$
\begin{bmatrix}
    \vec{y}_{m+1} & \vec{G} \hat{\vec{C}}
\end{bmatrix}
=
\vec{G}\vec{C}
\sim \operatorname{Gaussian}(d,d-m).
$$
Moreover, since $\vec{G}$ was independent of $\vec{X}$ and $\vec{Y}$, so is $\vec{G}\vec{C}$.
Thus, $\vec{G}\hat{\vec{C}}$ is independent of $\vec{X}$, $\vec{Y}$, and $\vec{x}_{m+1}$ and $\vec{y}_{m+1}$.

Now, using the orthogonality of $\vec{C}$,
$$
\vec{A} 
= 
\begin{bmatrix}\vec{Y} & \vec{G}\vec{C} \end{bmatrix} \begin{bmatrix}\vec{X}^\T \\ \vec{C}^\T\vec{Z}^\T \end{bmatrix}
=
\begin{bmatrix}
    \vec{Y} & \vec{y}_{m+1} & \vec{G}\hat{\vec{C}}
\end{bmatrix}
\begin{bmatrix}
\vec{X}^\T \\
\vec{x}^\T \\
(\vec{Z}\hat{\vec{C}})^\T
\end{bmatrix}
$$
Relabeling quantities gives the result.
:::

We will also use the following fact:


:::{prf:lemma}
Suppose $a_1 + \cdots + a_d = d/2$ and $a_i \leq 1$.
Then at least $d/4$ of the $a_i$ are larger than $1/4$.
:::

:::{prf:proof}
:class: dropdown
:enumerated: false
Let $L = \{ i : |a_i| > 1/4 \}$.
Suppose, for the sake of contradiction, that $|L| < d/4$.
Then
$$
\frac{d}{2} 
= \sum_{i\in L} a_i + \sum_{i\in L^{\textsf{c}}} a_i
\leq |L| \cdot 1 + (d-|L|) \cdot \frac{1}{4}
<  \frac{d}{4} \cdot 1 + d \cdot \frac{1}{4}
= \frac{d}{2}.
$$
This is a contradiction, so the result holds.
:::


:::{prf:proof} @thm:lower-diag
:enumerated: false

Set 
$$
d = \frac{1}{\varepsilon}
,\qquad
m = \frac{d}{4},
\qquad
k = d-m,
$$
and let $\vec{A} \sim\operatorname{Gaussian}(d)$.


First, note that $\|\vec{A}\|_\F^2$ is a Chi-squared random variable with $d^2$ degrees of freedom, and therefore has mean $d^2$.
Therefore, by Markov's inequality, 
$$
\PP\Big[ \|\vec{A}\|^2 \leq 50 d^2  \Big] \geq \frac{49}{50}.
$$

The algorithm want's to approximate the diagonal of $\vec{A} = \vec{Y}\vec{X}^\T + \vec{G}\vec{Z}^\T$, but after $m$ queries it only knows $\vec{X}$, $\vec{Y}$, and $\vec{Z}$.
We will show the diagonal of $\vec{G}\vec{Z}^\T$ is too random, so that it can't be approximated well (in squared norm relative to $O(d)$).

Write the rows of $\vec{G}$ and $\vec{g}_i^\T$ and the columns of $\vec{Z}^\T$ as $\vec{z}_i$. 
Then,

$$
\vec{d} = \operatorname{diag}(\vec{G}\vec{Z}^\T)
= 
\big[ \vec{g}_1^\T \vec{z}_1 , \vec{g}_2^\T\vec{z}_2, \cdots ,\vec{g}_d^\T \vec{z}_d \big].
$$
Note that since the entries of $\vec{G}$ are all independent Gaussians, the entries of $\vec{d}$ are also independent Gaussians! 
In particular, 
$$
[\vec{d}]_{\,i} 
= \vec{g}_i^\T \vec{z}_i 
= \sum_{j=1}^{k} [\vec{g}_i]_{\,j} \cdot [\vec{z}_i]_{\,j}
\sim \sum_{j=1}^{k} \mathcal{N}(0,[\vec{z}_i]_{\,j}^2)
\sim \mathcal{N}(0,\|\vec{z}_i\|_2^2).
$$
The algorithm knows this, so the best choice of a guess for $\vec{d}$ is $\vec{d} = \vec{0}$; indeed, each entry is a mean zero Gaussian.
Thus, the algorithm can output a guess $\widetilde{\vec{d}} = \operatorname{diag}(\vec{Y}\vec{X}^\T)$ for $\operatorname{diag}(\vec{A})$.
Then 
$$
\|\operatorname{diag}(\vec{A}) - \widetilde{\vec{d}} \|_2^2
= \|\vec{d}\|_2^2.
$$
We will show that $\|\vec{d}\|_2^2$ is typically large, so that the output always has high expected error.
Intuitively this makes sense; since the entries of $\vec{d}$ are Gaussian, and since the variance of a Gaussian is 1, 
$$
\mathbb{E}[\|\vec{d}\|_2^2]
= \sum_{i=1}^{d} \mathbb{E}[ |[\vec{d}]_{\,i}|^2 ]
= \sum_{i=1}^{d} \|\vec{z}_i\|_2^2
= \| \vec{Z} \|_\F^2
= k
= \frac{3}{4}d.
$$
This is on the right scale, since $\|\vec{A} - \vec{I}\circ\vec{A}\|_\F^2 \leq \|\vec{A}\|_\F^2 = O(d^2) = O(\varepsilon \:d)$.
However, an expectation bound doesn't prove a probability bound; it's possible to have a large expectation because of a very large and very unlikely value. 

In this particular case, we could simply use a [standard anti-concentration result for Gaussians](https://mathoverflow.net/a/95108).
However, for the sake of mirroring the proof we used in the paper, we will use an alternate approach.

Note that the $\|\vec{z}_i\|_2^2$ can't all be super small.
Indeed, let 
$$
L = \{ i : \|\vec{z}_i\|_2^2 \geq 1/4 \}.
$$
Since $[\vec{X}~\vec{Z}]$ is an orthogonal matrix $\|\vec{z}_i\|_2^2 \leq 1$ and $\|\vec{Z}\|_\F^2 = k \geq d/2$.
Thus, by Lemma 2, $|L| \geq d/4$.

Now, define 
$$
P = \{ i \in L : | [\vec{d}]_{\,i} | > 10^{-3} \}.
$$
If $i\in L$, just by standard properties of Gaussians,
$$
\mathbb{P}\big[ | [\vec{d}]_{\,i} | > 10^{-3} \big] \geq \frac{99}{100}.
$$
This implies that $P$ is not much smaller than $L$.
In particular, 
$$
\mathbb{P}[ |P| > |L|/2 ]
\geq \frac{49}{50}.
$$

We have now established that 

- $|P| > |L|/2 \geq d/8$ with probability at least $49/50$, and 
- $\|\vec{A}-\vec{I}\circ\vec{A}\|_\F^2 \leq \|\vec{A}\|_\F^2 \leq 50 d^2$ with probability at least $49/50$.

Thus, by a union bound, both events happen simultaneously with probability at least $24/25$. 

In this case, for some constant $D$, 
$$
\| \vec{I}\circ \vec{A} - \widetilde{\vec{A}} \|_\F^2
\geq
\| \vec{d} \|_2^2
= \sum_{i=1}^{d} [\vec{d}]_{\,i}^2
\geq \frac{d}{4} \cdot (10^{-3})^2
\geq D \varepsilon \: \| \vec{A} - \vec{I}\circ\vec{A} \|_\F^2.
$$
Replacing $\varepsilon$ with $\varepsilon/D$ throughout the proof gives the result with $C = 1/(4D)$.

:::


### What changes in the general case?

The above proof actually captures pretty much all of the important ideas we use for the lower-bound provided in the paper. 
There are a number of details which differ. 
The most important are perhaps the following:

- In the above proof, the entries of $\vec{d}$ were independent. For other sparsity patterns this won't be the case. 
While we used independence in our computation of the variance (which was just for intuition), we did not use independence in the actual proof.
- In the paper, we use the distribution $\vec{G}^\T\vec{G}$. 
This somewhat complicates the analog of Lemma 1, as well as computations (since things aren't all Gaussian anymore). 
However, it gives a lower-bound against the larger class of algorithms which can perform matvecs with the transpose.
- In the paper we are more careful about saying things like "let $d = 1/\varepsilon$" and "$\vec{d} = \vec{0}$ is the best choice", and handle these cases explicitly / precisely.