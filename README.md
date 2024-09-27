# NonlinearNormalForm

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://bmad-sim.github.io/NonlinearNormalForm.jl/stable/)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://bmad-sim.github.io/NonlinearNormalForm.jl/dev/)
[![Build Status](https://github.com/bmad-sim/NonlinearNormalForm.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/bmad-sim/NonlinearNormalForm.jl/actions/workflows/CI.yml?query=branch%3Amain)

## Overview

_This package is still in development, and certain bugs/feature and syntax changes must be expected._

This package provides routines for doing canonical perturbation theory on nonlinear Hamiltonian (and nearly Hamiltonian) maps, including parameters, using Lie algebraic methods. Functions are provided which, given some dynamical map expressed as a truncated power series in the variables (and parameters), can be used for calculating and analyzing the canonical transformation of the variables to the _**normal form**_ - coordinates where the nonlinear motion lies on action-dependent circles in phase space. This allows for easy calculation of invariants and other (parameter-dependent) properties of the map. Optionally, one single resonance may be left in the map, leaving a single resonance normal form from which resonance driving terms and the positions of fixed points in phase space may calculated. This package may be of particular interest to those in accelerator physics, celestial mechanics, electron microscopy, geometrical optics, and plasma physics.

## Setup

To use this package, in the Julia REPL run:

```julia
import Pkg; Pkg.add(url="https://github.com/bmad-sim/NonlinearNormalForm.jl")
```

## Basic Usage

This package imports and reexports [`GTPSA.jl`](https://github.com/bmad-sim/GTPSA.jl), a library for computing real and complex truncated power series to arbitrary orders in the variables and parameters. Before using `NonlinearNormalForm.jl`, you should have some familiarity with `GTPSA.jl`. 

The package currently provides various functionalities already provided by the [Full Polymorphic Package (FPP)](https://github.com/bmad-sim/bmad-ecosystem/blob/main/forest/fpp_manual/fpp-manual.pdf) written in Fortran90. This includes real and complex differential algebraic maps with properly overloaded operators, map composition and inversion (using routines provided by `GTPSA.jl`), parametric normal form calculation routines optionally including a "coasting" plane (including a constant "energy-like" canonical variable), factorization of the normalizing map, Lie operators including a quaternion for spin, calculations such as `exp` of Lie operators to construct Lie maps or the `log` of Lie maps to obtain the Lie operator, and one resonance normal form analysis tools.

After a `DAMap` is calculated via polymorphic tracking of the truncated power series, the map can be analyzed using the routines here. Some example maps randomly generated by FPP are provided in the `test` directory.

```julia
julia> m = read_fpp_map("test/spin_res/test.map") # read one of the randomly-generated maps from FPP into Julia

julia> m_lin = cutord(m, 2); # extract the linear part in orbital

julia> m_nonlinear = inv(m_lin) ∘ m; # remove the linear part

julia> F = log(m_nonlinear); # Get the Lie operator (including quaternion) generating nonlinear part

julia> m = m_lin ∘ exp(F);  # Reconstruct same map using Lie exponent and linear part separately

julia> a = normal(m);  # Calculate the nonlinear (parametric) normalizing canonical transformation

julia> R_z = inv(a) ∘ m ∘ a; # Nonlinear amplitude-dependent rotation in regular phase space (x, px, …)

julia> c = to_phasor(m);  # Get the transform to phasors basis √(J)*exp(±im*ϕ)

julia> R_J = inv(c) ∘ R_z ∘ c; # Nonlinear amplitude-dependent rotation in phasors basis

julia> a_spin, a0, a1, a2 = factorize(a); # Spin part, nonlinear parameter-dependent fixed point, a1, a2

julia> Σ = equilibrium_moments(m, a); # Calculate equilibrium sigma matrix when fluctuation-dissipation

julia> a = normal(m, m=[0; 1], m_spin=[-1]);  # Leaving in a Q_y - Q_spin resonance
```

## Acknowledgements
Thanks to Etienne Forest, the creator of FPP, for his significant time and patience in teaching the normal form methods and guiding the implementation of the tools in this package.
