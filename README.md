# RationalPolygons.jl

[![](https://img.shields.io/badge/docs-stable-blue.svg)](https://justus-springer.github.io/CStarSurfaces.jl/stable)
[![](https://img.shields.io/badge/docs-dev-blue.svg)](https://justus-springer.github.io/CStarSurfaces.jl/dev)

RationalPolygons.jl is a pure Julia package for computations with rational convex polygons. It implements

- counting lattice points, Ehrhart Theory, normal forms, automorphism groups,
  computation of subpolygons,
- various classification algorithms] for integral and rational polygons.

I have written RationalPolygons.jl in the span of about six months while
working on a [joint
project](https://justus-springer.github.io/RationalPolygons.jl/dev/#BS24) with
Martin Bohnert. Its main purpose is to provide reference implementations of the
classification algorithms developed in our paper. However, it also implements
many more basic algorithms for computations with rational polygons that I
believe might be useful in other projects. RationalPolygons.jl does not make
use of any external computer algebra system but implements all necessary
algorithms, including two-dimensional euclidian geometry, from scratch in pure
Julia. This allows for quite good performance, with computations involving
billions of polygons being feasable on a personal computer.

# Documentation

[The full documentation can be found here.](https://justus-springer.github.io/CStarSurfaces.jl/dev)
