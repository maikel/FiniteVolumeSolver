[![Build on Ubuntu 22.04](https://github.com/maikel/FiniteVolumeSolver/actions/workflows/ubuntu.yml/badge.svg)](https://github.com/maikel/FiniteVolumeSolver/actions/workflows/ubuntu.yml)

# Finite Volume Solver

This C++17 project provides a framework to solve hyperbolic equations with finite volume methods.
The framework assumes an external structured grid library like AMReX or SAMRAI to adaptively refine regions of interest.

One of the main design philosophies is to make the numerical methods re-usable and consistent across different grid implementations.
We provide some standard flux methods which can be used out-of-the-box with a large set of hyperbolic equations.
Furthermore, a lot of helper classes exist to generate distributed grids for data management which simplifies and modernizes the usage of libraries like AMReX or SAMRAI.

At last, this library is also capable of handling embedded boundaries in a dimensionally split setting as defined in [Klein2009].

We provide a [user documentation](http://page.mi.fu-berlin.de/ghastermann/fvs-agklein-doc/) and a [Doxygen documentation](http://page.mi.fu-berlin.de/ghastermann/fvs-agklein-dox/).
