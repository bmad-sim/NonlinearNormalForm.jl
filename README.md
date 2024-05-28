# NonlinearNormalForm

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://bmad-sim.github.io/NonlinearNormalForm.jl/stable/)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://bmad-sim.github.io/NonlinearNormalForm.jl/dev/)
[![Build Status](https://github.com/bmad-sim/NonlinearNormalForm.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/bmad-sim/NonlinearNormalForm.jl/actions/workflows/CI.yml?query=branch%3Amain)

## Overview

_This package is still in an early development stage_

This package provides routines for calculating parameter-dependent normal forms of nonlinear Hamiltonian (and nearly Hamiltonian) maps expressed as truncated power series in the variables, using Lie algebraic methods. Basically, given some map, it computes the canonical transformation of the variables to coordinates where the nonlinear motion lies on circles in phase space. This allows calculation of invariants of the motion, and when leaving one resonance in the map, a single-resonance normal form. This package may be of particular interest to those in accelerator physics, electron microscopy, and astronomy.

While transfer maps act on phase space variables, the useful Lie representation of the map acts on functions of the phase space variables; e.g. for some Hamiltonian map $\vec{\zeta}(\vec{x})$, the Lie map $\mathcal{M}\_{\vec{\zeta}}$ acts as $$\mathcal{M}\_{\vec{\zeta}}f=f\circ\vec{\zeta} \ .$$  

$\mathcal{M}_{\vec{\zeta}}$ is a [Koopman operator](https://en.wikipedia.org/wiki/Composition_operator). This replaces the nonlinear problem of a finite dimension (describing how the phase space variables transform to infinite order in the phase space variables), to a linear problem of infinite dimension (describing how each monomial coefficient in the Hamiltonian map transforms). The number of monomial coefficients is simply truncated at some order, chosen by the user.

Because the problem is now linear, the transformation of the monomial coefficients could be expressed using a (truncated) matrix, however this is notoriously slow and inconvenient to work with. A Lie algebraic formulation is much more elegant, simpler, and faster to use.

Hamilton's equations expressed using the Poisson bracket are

$$ \frac{d}{dt} q = -\lbrace H, q\rbrace\ , \ \ \ \ \ \  \frac{d}{dt} p = -\lbrace H, p\rbrace \ . $$

For simplicity, we define an operator using the Poisson bracket as 
$$:f:g \triangleq \lbrace f, q \rbrace$$

The Lie map and Hamiltonian map for an infinitesimal timestep $dt$ are respectively expressed as 

$$\mathcal{M}\_{t_0\rightarrow t_0+dt} = \exp{(:-dt\ H:)} \ , \ \ \ \ \ \ \vec{\zeta}\_{t_0\rightarrow t_0+\Delta t}= \exp{(:-dt\ H:)}\vec{x}$$

where $H$ is the Hamiltonian, a scalar function. While using the Hamiltonian as a generator of the time evolution (canonical transformation) in this way is useful, it is not always practical. In a particle tracking code for example, the motion may not always be exactly symplectic, e.g. when there is some small dissipation due to synchrotron radiation emission. In this case, a Hamiltonian cannot be obtained from an integrator, however the force field $\vec{F}$ on the particle is of course known from the integrator. Using the operator isomorphism $\vec{F}\cdot\vec{\nabla} \triangleq \ :-H: $, we can instead write

<!--- $$\mathcal{M}\_{t_0\rightarrow t_0+dt} = \exp{(\vec{F}\cdot\vec{\nabla}\ dt)} \ , \ \ \ \ \ \ \vec{\zeta}\_{t_0\rightarrow t_0+\Delta t}= \exp{(:-H\ dt:)}\vec{x}$$

An isomorphism can be used to express the above equations, which use the Hamiltonian, to instead use the force field




a _Lie operator_ generating the map. If the system is Hamiltonian, the generator of the map would simply be a Poisson bracket, however vector fields are used so the formalism applies correctly even for nearly-Hamiltonian maps. --->


## Setup

To use this package, in the Julia REPL run:

```juliai
import Pkg; Pkg.add(path="https://github.com/bmad-sim/NonlinearNormalForm.jl")
```

## Basic Usage

This package imports and reexports [`GTPSA.jl`](https://github.com/bmad-sim/GTPSA.jl), a library for computing real and complex truncated power series to arbitrary orders in the variables and parameters. Before using `NonlinearNormalForm.jl`, you should have some familiarity with `GTPSA.jl`. 


