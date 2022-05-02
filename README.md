# PathIntersections

[![Build status](https://github.com/cgt3/PathIntersections.jl/workflows/CI/badge.svg)](https://github.com/cgt3/PathIntersections.jl/actions)
[![Build Status](https://travis-ci.com/cgt3/PathIntersections.jl.svg?branch=main)](https://app.travis-ci.com/github/cgt3/PathIntersections.jl)
[![Build Status](https://ci.appveyor.com/api/projects/status/github/cgt3/PathIntersections.jl?svg=true)](https://ci.appveyor.com/project/cgt3/PathIntersections-jl)
[![Coverage](https://codecov.io/gh/cgt3/PathIntersections.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/cgt3/PathIntersections.jl)
[![Coverage](https://coveralls.io/repos/github/cgt3/PathIntersections.jl/badge.svg?branch=main)](https://coveralls.io/github/cgt3/PathIntersections.jl?branch=main)


## Introduction
This package locates intersections between parameterized curves (paths) and Cartesian meshes. The package locates intersections by walking along the curve using the user provided step size, `ds`, and sensing for boundary crossings. When a boundary is crossed, a hybrid bracketed secant-bisection method is used to determine the location of the intersection within the user-specified tolerances `arc_tol` and `single_tol`. 

While the use of a bracketed method is robust in many situations, it does not allow all non-simple intersections to be found; cases that are currently unsupported or may result in unsatisfactory accuracy are noted below.

## Usage and Examples

### Single Curve

### Multiple Curves

### Tolerance Defitions
- `arc_tol`
The main tolerance used to find intersections. This tolerance is applied to the Euclidean distance between the upper and lower bounds of the secant-bisection method at each iteration. The method stops when the Euclidean distance between these points is less than or equal to `arc_tol`.

<img src="/figures/tolerance_arcTol.png" width="500">


- `single_tol`
The secondary tolerance for corner intersections. This tolerance is used to determine if another dimension participated in previously found intersection (i.e. to determine if the intersection went through a corner). It is applied to the distance between the previously estimated intersection point and the nearest mesh edge in the dimension being checked. It is, with the exception of folded intersections, a weaker condition than `arc_tol`.

<img src="/figures/tolerances_corners.png" width="500">

### Assumptions
1. Each array in the array of mesh coordinates, `coords`,
    1. is sorted in increasing order (i.e. `coords[dim][i] < coords[dim][i+1]`),
    2. does not contain duplicates, and
    3. contains the start and end points of the domain in that dimension.
2.  All curves start at `s = 0` and end at `s = 1`. These points will be checked regardless of the step size provided.
3.  The provided step size is sufficiently small in comparison to the mesh for intersections to be sensed.

### Notes:
1. This package is only concerned with mesh-curve intersections, not curve-curve intersections.
2. Intersections are returned sorted in increasing order according to their `s` values.
3. Curves may start and arbitrarily enter and leave regions outside the domain defined by the mesh.

## Supported Cases
### Simple, Single-Dimension Intersections

<img src="/figures/intersection_simpleSingleDim.png" width="500">

### Simple, Multiple-Dimension Intersections

<img src="/figures/intersection_simpleMultiDim.png" width="500">

### Corners

<img src="/figures/intersection_corner.png" width="500">

## Unsupported and Degenerate Cases
### Unsupported: General Non-simple Intersections

<img src="/figures/intersections_nonSimple.png" width="500">

### Unsupported: Bad Step Size

<img src="/figures/intersections_badStepSize.png" width="500">

### Degenerate Input: Folded Intersections

<img src="/figures/tolerance_foldedIntersection.png" width="500">
