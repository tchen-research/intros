---
title: Randomized Algorithms
description: What is a randomized algorithm?
---


We would like the algorithms we use to produce exactly the right answer all of the time. 
However, in many cases, this is not practically possible, as an algorithm would have to brute-force through all possible solutions to find the right one.
Randomized algorithms are a class of algorithms that use randomness as a tool to help find good (but not necessarily exact) solutions quickly.

## Toy Example

Here's a toy example: Suppose we want to find the best restaurant in NYC.
To do this, we'd have to try every single restaurant in NYC, an impossible task.
However, if instead we're just asked to find a restaurant in the top 10% of all restaurants, we could simply try a random selection of a small number of restaurants (e.g. 20), and report the best one we found.

Let's try to analyze this problem mathematically. Suppose there are $N$ restaurants (in real life $N\approx 20000$). 
Each restaurant has a quality score which is unknown, until we visit it (for simplicity, we'll assume that these scores are unique and deterministic).

How can we find a restaurant in the top 10% of all restaurants?

Here's a natural algorithm: 
We'll commit to visiting $M$ restaurants that we will pick completely randomly (i.e. by putting the names of all restaurants in a giant hat and drawing $M$ of them).
Then, we will visit each of these $M$ restaurants, and return the one with the highest score.
The big question is: *how large does $M$ have to be to guarantee that we will find a restaurant in the top 10% of all restaurants?*


Unless $M$ is more than 90% of $N$, then there is a chance that we will not find a restaurant in the top 10%; we might get exceptionally unlucky and draw only restaurants in the bottom 90%. 
However, because we've chosen the restaurants randomly, we know this is very unlikely! 

Let's call the set of restaurants in the bottom 90% of all restaurants the "bad restaurants". And to simplify analysis, let's assume that when we draw a restaurant from the hat, we return it to the hat before drawing the next one (so we might draw the same restaurant multiple times).

- Each time we draw a restaurant, the probability that we draw a bad restaurant is $0.9$.
- The chance we only draw bad restaurants after $M$ draws is $0.9^M$.
- The chance we draw at least one good restaurant after $M$ draws is $1 - 0.9^M$.

We can now choose $M$ to get the desired probability of success. 
For instance, suppose we want this procedure to succeed with $95\%$ probability.
Then we can set $1 - 0.9^M \geq 0.95$ (or equivalently $0.9^M \leq 0.05$). 
Solving this equation (by taking the logarithm of both sides) we see that we need $M\geq 29$ to guarantee that we will find a restaurant in the top 10% of all restaurants with probability at least 95%.

### The general case

What if we want to find a restaurant in the top $\varepsilon$ fraction of all restaurants with probability at least $1-\delta$? 
Then we can set $1 - (1-\varepsilon)^M \geq 1-\delta$. Solving for $M$ we see that we need 
$$
M = \frac{\log(\delta)}{\log(1-\varepsilon)}
\approx \frac{\log(1/\delta)}{\varepsilon},
$$
where we have used that $\log(1-\varepsilon) \approx -\varepsilon$ for small $\varepsilon$ and $\log(\delta) = -\log(1/\delta)$.

So we see how the cost of the algorithm depends on the desired accuracy $\varepsilon$ and the desired probability of success $1-\delta$.
Interestingly, the dependence on the probability of success is very mild. 
For instance, if in our example we wanted to succeed at least 99.9% of the time, we would only need to increase $M$ from 29 to 66.
On the other hand, the dependence on the accuracy $\varepsilon$ is much worse. Finding a restaurant in the top 1% of all restaurants requires visiting roughly 10 times as many restaurants as finding a restaurant in the top 10% of all restaurants.

## The role of randomness

A deterministic algorithm is one does not rely on randomness at all. 
For instance, we might order all of the restaurants in alphabetical order, and then visit the first $M$ of these restaurants.
However, it is very hard to guarantee such an approach will work. 
For instance, what if all of the good restaurants happen to be later in the alphabet? 
While our intuition tells us this is unlikely, formalizing this intuition is very difficult even in this simple example (and may not even be possible in other more complicated examples).

On the other hand, the use of randomness allowed us to provide strong guarantees about the performance of our very simple algorithm.

