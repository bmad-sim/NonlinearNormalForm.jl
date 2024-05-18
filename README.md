# NonlinearNormalForm

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://bmad-sim.github.io/NonlinearNormalForm.jl/stable/)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://bmad-sim.github.io/NonlinearNormalForm.jl/dev/)
[![Build Status](https://github.com/bmad-sim/NonlinearNormalForm.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/bmad-sim/NonlinearNormalForm.jl/actions/workflows/CI.yml?query=branch%3Amain)

## Overview

This package provides routines for calculating parameter-dependent normal forms of nonlinear Hamiltonian (and nearly Hamiltonian) maps expressed as truncated power series in the variables, using Lie algebraic methods. Basically, given some map, it computes the canonical transformation of the variables to coordinates where the nonlinear motion lies on circles in phase space. 

While transfer maps act on phase space variables, the useful Lie representation of the map acts on functions of the phase space variables; e.g. for some Hamiltonian map $\vec{\zeta}(\vec{x})$, the Lie map $\mathcal{M}\_{\vec{\zeta}}$ acts as $$\mathcal{M}\_{\vec{\zeta}}f=f\circ\vec{\zeta} \ .$$  

$\mathcal{M}_{\vec{\zeta}}$ is a [Koopman operator](https://en.wikipedia.org/wiki/Composition_operator). This replaces the nonlinear problem of a finite dimension (describing how the phase space variables transform to infinite order in the phase space variables), to a linear problem of infinite dimension (describing how each monomial coefficient in the Hamiltonian map transforms). The number of monomial coefficients is simply truncated at some order, chosen by the user.

Because the problem is now linear, the transformation of the monomial coefficients could be expressed using a (truncated) matrix, however this is notoriously slow and inconvenient to work with. A Lie algebraic formulation is much more elegant, simpler, and faster to use.

We define 
Hamilton's equations expressed using the Poisson bracket are

$$ \frac{d}{dt} q = -\lbrace H, q\rbrace\ , \ \ \ \ \ \  \frac{d}{dt} p = -\lbrace H, p\rbrace \ . $$

For simplicity, we now define an operator using the Poisson bracket: $:f:g \triangleq \lbrace f, q \rbrace$

The Lie map and Hamiltonian map are respectively expressed as 

$$\mathcal{M} = \exp{(\hat{F})} \ , \ \ \ \ \ \ \vec{\zeta}= \exp{(\hat{F})}\vec{x}$$

where $\hat{F}$ is a _Lie operator_ generating the map. If the system is Hamiltonian, the generator of the map would simply be a Poisson bracket, however vector fields are used so the formalism applies correctly even for nearly-Hamiltonian maps.


## Setup

To use this package, in the Julia REPL run:

```julia
import Pkg; Pkg.add(path="https://github.com/bmad-sim/NonlinearNormalForm.jl")
```
