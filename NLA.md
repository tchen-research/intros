---
title: Numerical Linear Algebra
description: Brief introduction to numerical linear algebra
---


Numerical Linear Algebra is about algorithms to "do linear algebra on computers".

If you haven't seen linear algebra in a while, it's just a certain collection of types of math problems that includes things like [linear systems of equations](https://en.wikipedia.org/wiki/Linear_system) such as $3x+2y^2=7$ and $2x-y= 0$.
Linear algebra has a ton of applications, in part because many complex algorithms rely on solving linear algebra problems as subroutines.

The above system only has two variables, so you can probably figure out the solution by just playing around.
However, once more variables are involved, we need to use some kind of well-defined procedure; i.e. an algorithm!
For example, you might have learned about things like ["row-reduction"](https://en.wikipedia.org/wiki/Gaussian_elimination#Row_operations) or ["Gaussian elimination"](https://en.wikipedia.org/wiki/Gaussian_elimination), which are algorithms (that you can run by hand) to solve linear systems in a principled way.
Since doing calculations is time-consuming, we would like to come up with algorithms that give us the solution in as few operations as possible.

Of course, even with the best algorithms, we might run into problems too large to solve by hand. 
In fact, many modern applications require solving systems with *millions* of variables, which we definitely can't solve by hand (imagine how long it would take to even write down a million numbers, let alone do arithmetic with them).
Fortunately, computers are really good at doing lots of basic tasks quickly!

However, once we start trying to do math on computers, we run into some difficulties. 
Computers are discrete and finite; they work with "1"s and "0"s and clearly have a limited amount of storage/processing power/etc.
On the other hand, there are an infinite number of numbers! 
To get around this issue, it is common to use an approximate number system (e.g. floating point numbers) which can approximately represent a large range of numbers. 

While this enables us do math, it has the undesirable consequence that every time we do basic arithmetic operations, such as adding two numbers, we might make a tiny error.
After all the true answer might not even be representable in our number system.
Fortunately, the number systems that are commonly used are good enough that we know that every error we make is tiny.
However, while each basic operation only results in a tiny error, we need to make sure that we design algorithms in a way that avoids all these tiny errors compounding into a large error (or at least make sure we know when this happens).


