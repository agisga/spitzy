---
layout: post
title:  "Documentation"
date:   2015-05-05 13:10:08
---

# Table of Contents

1. [Introduction](#introduction)

2. [Initial Value Problems](#ivp)

  * [Forward Euler](#euler)

  * [Adams-Bashforth of Order 2](#ab2)

  * [Dormand-Prince](#dopri)

3. [Two-point Boundary Value Problems](#bvp)

  * [Linear Finite Element Galerkin](#lin_fin_elt)

4. [2D Poisson's Equation](#poissons_eq)

  * [Five-point Laplacian](#five_pt)

5. [1D Linear Advection Equation](#advection_eq)

  * [Upwind](#upwind)

  * [Leapfrog](#leapfrog)

  * [Lax-Friedrichs](#lax_friedrichs)

  * [Lax-Wendroff](#lax_wendroff)

6. [Conclusions and Future Work](#conclusion)

7. [References](#references)

---------------------------------------------------

<div id='introduction'/>
##Introduction
TO BE WRITTEN

<div id='ivp'/>
##Initial Value Problems
TO BE WRITTEN

<div id='euler'/>
###Forward Euler
TO BE WRITTEN

<div id='ab2'/>
###Adams-Bashforth of Order 2
TO BE WRITTEN

<div id='dopri'/>
###Dormand-Prince
TO BE WRITTEN

<div id='bvp'/>
##Two-point Boundary Value Problems

A two-point boundary value problem that can be solved with 'spitzy' is an ordinary differential equation of the following general form:

* ODE: $-(\alpha u')'(x) + (\beta u')(x) + (\gamma u)(x) = f(x)$, for $a < x < b$,
* where $\alpha$, $\beta$ and $\gamma$ are constants or continuous functions of $x$ on $[a,b]$,
* with Dirichlet boundary conditions: $u(a) = u\subscript{a}$ and $u(b) = u\subscript{b}$.

Currently the linear finite element Galerkin method is the only scheme implemented in `spitzy` to solve this problem.

<div id='lin_fin_elt'/>
### Linear Finite Element Galerkin

This section is largely based on section 12.4 in [QSS07]. The interested reader is referred there for a more rigorous presentation of what follows.

#### Background

Consider an ordinary differential equation of the following general form:

* ODE: $-(\alpha u')'(x) + (\beta u')(x) + (\gamma u)(x) = f(x)$, for $a < x < b$,
* where $\alpha$, $\beta$ and $\gamma$ are continuous functions of $x$ on $[a,b]$,
* with Dirichlet boundary conditions: $u(a) = u\subscript{a}$ and $u(b) = u\subscript{b}$.

We are looking for a [weak solution](http://en.wikipedia.org/wiki/Weak_solution) to this problem, that is, a function $u$ that satisfies the integral equation
$$\underset{[a,b]}{\int} \alpha u' v' dx + \underset{[a,b]}{\int} \beta u' v dx + \underset{[a,b]}{\int} \gamma u v dx = \underset{[a,b]}{\int} f v dx,$$
for every $v$ in the so-called [test-function space](http://en.wikipedia.org/wiki/Distribution_%28mathematics%29#Test_functions_and_distributions) $V$, which basically contains functions that have square-integrable [distributional derivatives](http://en.wikipedia.org/wiki/Distribution_%28mathematics%29#Derivatives_of_distributions) and vanish on the boundary of the domain of the ODE. If $u$ is a solution to the original formulation of the ODE then it also satisfies to integral formulation. However, a solution $u$ of the integral equation might not be twice differentiable.

#### Galerkin Method

The [Galerkin method](http://en.wikipedia.org/wiki/Galerkin_method) approximates $V$ with a finite dimensional functional space $V\subscript{h}$, which leads to a finite number of test functions $v$ that need to be tested against $u$ with the above integral equation. We denote the basis functions of $V\subscript{h}$ by $\phi\subscript{1}, \phi\subscript{2}, \ldots, \phi\subscript{N}$. Moreover, if we assume that $u$ also lies in the space $V\subscript{h}$, then the problem reduces to a linear system
$$A\vec{u} = \vec{f},$$
where $A$ has entries
$$a\subscript{ij} = \underset{[a,b]}{\int} \alpha \phi\subscript{i}' \phi\subscript{j}' dx + \underset{[a,b]}{\int} \beta \phi\subscript{i}' \phi\subscript{j} dx + \underset{[a,b]}{\int} \gamma \phi\subscript{i} \phi\subscript{j} dx,$$
and the right hand side vector $\vec{f}$ has entries
$$f\subscript{i} = \underset{[a,b]}{\int} \alpha f' \phi\subscript{i}' dx + \underset{[a,b]}{\int} \beta f' \phi\subscript{i} dx + \underset{[a,b]}{\int} \gamma f \phi\subscript{i} dx.$$

The structure of $A$ and the degree of accuracy of the numerical solution depends on the form of the basis functions $\phi\subscript{1}, \phi\subscript{2}, \ldots, \phi\subscript{N}$, that is, on the choice of $V\subscript{h}$.

#### Finite Element Method

The [finite element method](http://en.wikipedia.org/wiki/Finite_element_method) chooses $V_h$ to be the space of continuous piecewise polynomials which are defined on subintervals of $[a, b]$ and vanish at $a$ and $b$.

Here, we only consider polynomials of degree 1, that is, continuous piecewise linear functions.
Given an equally spaced grid $a = x\subscript{0} < x\subscript{1} < \ldots < x\subscript{n-1} < x\subscript{n} = b$, define for $i=1,2,\ldots,n-1$ the function $\phi\subscript{i}$ to be the continuous piecewise linear function that is equal to 1 at the node $x\subscript{i}$ and 0 at all other nodes. Thus, the basis functions have the following shape (image from a [wikipedia article](http://en.wikipedia.org/wiki/Finite_element_method), the red line represents a linear combination of the basis functions):

![Piecewise linear basis functions](/images/fin_elt_basis_func.svg?raw=true "Piecewise linear basis functions")

#### Implementation

The linear finite element Galerkin method is implemented in class `Bvp` of the Ruby gem [spitzy](https://github.com/agisga/spitzy). It takes as inputs the interval of $x$ values as well as the Dirichlet boundary conditions at the two edges of the interval, the desired number of equally spaced grid points on which the numerical solution will be evaluated, and the functions $\alpha(x)$, $\beta(x)$, $\gamma(x)$ and $f(x)$ as `Proc` objects (or as `Numeric` if the function is constant).

The linear finite element Galerkin method boils down to a linear system
$A\vec{u} = \vec{f}$, where the matrix $A$ is tridiagonal. Currently, the Ruby implementation in `spitzy` uses the `#solve` method from the `NMatrix` gem.
However, the memory usage as well as the speed of the algorithm can be improved, as soon as pull request [#301](https://github.com/SciRuby/nmatrix/pull/301), which implements a `#solve_tridiagonal` method, is merged into the project.

#### Examples

Now let's look at some examples.

<!-- * ODE: $-(\alpha u')'(x) + (\beta u')(x) + (\gamma u)(x) = f(x)$, for $a < x < b$,
* where $\alpha$, $\beta$ and $\gamma$ are continuous functions of $x$ on $[a,b]$,
* with Dirichlet boundary condition: $u(a) = u\subscript{a}$ and $u(b) = u\subscript{b}$. -->

#### Example 1. Constant coefficients

First we look at example, where $\alpha$, $\beta$ and $\gamma$ are constants independent of $x$:

* ODE: $-800\pi u'' + 8\pi u = 0$, $0 < x < 100$,
* BC: $u(0) = 10$, $u(100) = \frac{10}{\cosh(10)}$.

We compute the numerical solution:

```Ruby
require 'spitzy'
bvp_sol = Bvp.new(xrange: [0.0, 100.0], mx: 100, bc: [10.0, 0.00090799859], 
                  a: 800.0*Math::PI, b: 0.0, c: 8.0*Math::PI, f: 0.0)
```

The exact solution is 
$$u(x) = 10 \frac{\cosh\left(\frac{100 - x}{10}\right)}{\cosh(10)}.$$

We compare the obtained numerical to the exact solution by plotting both using the Ruby `gnuplot` gem:

![BVP example 1 plot](/images/bvp1.png?raw=true "BVP example 1 plot")

#### Example 2. Non-constant coefficients

Now, let's consider an example where $\alpha$, $\beta$, $\gamma$ and $f$ are all functions of $x$.

* ODE: $-\frac{d}{dx}(-\cos(x)u'(x)) + \sin(x)u'(x) + \cos(x) u(x) = -2\cos^2(x)$, $0 < x < 10$,
* BC: $u(0) = 0$, $u(10) = -9\sin(10)$.

We use `spitzy` to compute the numerical solution:

```Ruby
a = Proc.new { |x| -Math::cos(x) }
b = Proc.new { |x| Math::sin(x) }
c = Proc.new { |x| Math::cos(x) }
f = Proc.new { |x| -2.0*(Math::cos(x))**2 }
xrange = [0.0, 10.0]
bc = [0.0, (1.0 - 10.0) * Math::sin(10.0)]
bvp_sol = Bvp.new(xrange: xrange, mx: 100, bc: bc, 
                  a: a, b: b, c: c, f: f)
```

By taking derivatives one can check that the exact solution is given by
$$u(x) = (1-x)\sin(x).$$

Again, we plot both, the numerical and the exact solution, and observe that they approximately agree.

![BVP example 2 plot](/images/bvp2.png?raw=true "BVP example 2 plot")

<div id='poissons_eq'/>
##2D Poisson's Equation

Poisson's equation is the elliptic boundary value problem:

* $\Delta u = f$ in $\Omega$,
* with the boundary condition: $u=g$ on $\partial \Omega$.

Here we only consider the 2-dimensional version of Poisson's problem, where $\Delta u = \frac{d^2 u}{dx^2} + \frac{d^2 u}{dy^2}$.
We can discretise Poisson's equation by defining a mesh over the domain $\Omega$. If we denote by $\Omega\subscript{h}$ the set of the interior mesh points and by $\Gamma\subscript{h}$ the set of boundary mesh points, then the numerical solution needs to satisfy the discrete version of Poisson's problem written as:

* $\Delta\subscript{h} u = f$ in $\Omega\subscript{h}$,
* with the Dirichlet boundary condition: $u=g$ on $\Gamma\subscript{h}$.

Here $\Delta\subscript{h}$ is a consistent approximation of the operator $\Delta$, and it can be defined in various ways.

The only method currently implemented in `spitzy` defines the operator $\Delta\subscript{h}$ based on a five-point stencil, consisting of the point itself and its four immediate neighbors. We call it the five-point Laplacian method.

For more detail on the above and an introduction to various numerical methods for the solution of Poisson's equation we refer to [Arn11].

<div id='five_pt'/>
### Five-point Laplacian

In the current implementation of the five-point Laplacian method we assume a rectangular domain, i.e. $\Omega = [a,b]\times[c,d]$, and that each side of the domain is subdivided into subintervals of the same length $h$.

Under such conditions the following approximation to the Laplacian operator can be used
$$\Delta\subscript{h} u\subscript{i,j} = \frac{1}{h^2} (4u\subscript{i,j} - u\subscript{i+1,j} - u\subscript{i-1,j} - u\subscript{i,j+1} - u\subscript{i,j-1}),$$
where $x\subscript{i,j}$ denotes the point in the discretization $\Omega\subscript{h}$ corresponding to the $i$th coordinate in $x$-direction and the $j$th coordinate in $y$-direction, and $u\subscript{i,j}$ is an approximation to the true solution $u(x\subscript{i,j})$. After a moment of reflection it is clear that this amounts to the [second-order centered discretization of the second derivative](http://en.wikipedia.org/wiki/Finite_difference#Higher-order_differences) in both, the $x$- and the $y$-directions.

It is an easy exercise to see that this problem reduces to the linear system
$$A\vec{u} = \vec{b},$$
where $\vec{u} = (u\subscript{1,1}, u\subscript{2,1}, \ldots, u\subscript{m,1}, \ldots, u\subscript{1,k}, u\subscript{2,k}, \dots, u\subscript{m,k}, \dots, u\subscript{1,n}, u\subscript{2,n}, \dots, u\subscript{m,n})^T$, 
the matrix $A$ is pentadiagonal with constant values along each (sub-/super-)diagonal, and the right hand side $\vec{b}$ consists of values $f(x\subscript{i,j})$ with the boundary condition $g(x\subscript{i,j})$ subtracted when necessary. For more detail on the construction of $\vec{u}$, $A$ and $\vec{b}$ we refer to the definition of `Bvp` in the `spitzy` code.

#### Implementation

The five-point Laplacian method is implemented in class `Poissons_eq` of the Ruby gem [spitzy](https://github.com/agisga/spitzy). It takes as inputs the rectangular domain in form of a range in $x$-direction and a range in $y$-directions, the step size $h$, the Dirichlet boundary condition as a function of $x$ and $y$, and the right hand side function $f(x,y)$ (each function is passed as a `Proc` object, or as `Numeric` if the function is constant).

A mesh grid is generated with the method `#meshgrid` from the `NMatrix` gem. Currently this method is available only in the developement version of `NMatrix`.

Currently, the Ruby implementation in `spitzy` uses the `#solve` method from the `NMatrix` gem to solve the linear system.
Since the matrix $A$ is pentadiagonal other methods can significantly improve the computational speed as well as the memory usage. It would be straight forward to implement a modified [successive over-relaxation](http://en.wikipedia.org/wiki/Successive_over-relaxation) method, which would only perform multiplications corresponding to the non-zero elements of $A$ (it would even be unnecessary to save the matrix at all, because all entries of $A$ are either $-4/h^2$, $1/h^2$ or $0$). However, such a method should probably be implemented in C, because it is questionable if it would beat `NMatrix`'s `#solve` method (which is written in C) otherwise.

#### Examples

<div id='advection_eq'/>
##1D Linear Advection Equation
TO BE WRITTEN

<div id='upwind'/>
### Upwind
TO BE WRITTEN
 
<div id='leapfrog'/>
### Leapfrog
TO BE WRITTEN

<div id='lax_friedrichs'/>
### Lax-Friedrichs
TO BE WRITTEN

<div id='lax_wendroff'/>
### Lax-Wendroff
TO BE WRITTEN

<div id='conclusion'/>
##Conclusions and Future Work
TO BE WRITTEN

<div id='references'/>
##References

- [QSS07] A. Quarteroni, R. Sacco, F. Saleri (2007) *Numerical Mathematics*, 2nd ed., Texts in Applied Mathematics. Springer.
- [Arn11] Douglas N. Arnold (2011) *Lecture notes on Numerical Analysis of Partial Differential Equations*, version of 2011-09-05. Lecture notes MATH 8445 Numerical Analysis of Differential Equations, University of Minnesota. <http://www.ima.umn.edu/~arnold//8445.f11/notes.pdf>
