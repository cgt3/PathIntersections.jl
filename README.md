# PathIntersections

[![Build status](https://github.com/cgt3/PathIntersections.jl/workflows/CI/badge.svg)](https://github.com/cgt3/PathIntersections.jl/actions)
[![Build Status](https://travis-ci.com/cgt3/PathIntersections.jl.svg?branch=main)](https://app.travis-ci.com/github/cgt3/PathIntersections.jl)
[![Build Status](https://ci.appveyor.com/api/projects/status/github/cgt3/PathIntersections.jl?svg=true)](https://ci.appveyor.com/project/cgt3/PathIntersections-jl)
[![Coverage](https://codecov.io/gh/cgt3/PathIntersections.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/cgt3/PathIntersections.jl)
[![Coverage](https://coveralls.io/repos/github/cgt3/PathIntersections.jl/badge.svg?branch=main)](https://coveralls.io/github/cgt3/PathIntersections.jl?branch=main)


## Introduction
This package locates intersections between parameterized curves (paths) and Cartesian meshes in arbitrary dimension. The package locates intersections by walking along the curve using the user provided step size, `ds`, and sensing for boundary crossings. When a boundary is crossed, a hybrid bracketed secant-bisection method is used to determine the location of the intersection within the user-specified tolerances `arc_tol` and `single_tol`. 

While the use of a bracketed method is robust in many situations, it does not allow all non-simple intersections to be found; cases that are currently unsupported or may result in unsatisfactory accuracy are noted below.

## Usage and Examples
The package provides two functions, `find_path_intersections_single` and `find_path_intersections`, for finding intersections wrt a single curve or a set of curves, respectively. Both functions require the same arguments, differing only in their dimension. 

The arguments are:
- `coords`: A Cartesian mesh, which is to be provided in the form of an array of arrays with indexing `coords[dimension][edge index]`. This argument is the same for both `find_path_intersections_single` and `find_path_intersections` (i.e. both function only allow one mesh to be provided). Note the mesh does not need to be uniform.

- `curve`/`curves` : The curve/curves to be check against the mesh. The package expects curves to be callable using `pt = curve(s, curve_params...)` where `pt` is indexed `pt[dimension]`. For `find_path_intersections`, an array of curves with the indexing/signature: `pt = curves[curve index](s, curve_params...)`.

- `ds` : the/an array of step size(s) use for walking along the curve. The array is indexed `ds[curve index]`.
- `arc_tol` : the/an array of arc length tolerance(s) used as the stopping criteria for finding intersections. The array is indexed `arc_tol[curve index]`.
- `single_tol`: the/an array of single-dimension tolerance(s) used for identifying if other dimensions participated in an intersecton (i.e. whether the intersection occurred at a corner). The array is indexed `single_tol[curve index]`.
- `curve_params` : the/an array of parameters to be passed to the curve(s) along with `s`. While this is an optional argument for `find_path_intersections_single`, it is not optional for `find_path_intersections`, where an array indexed `curve_params[curve_index]` must be passed, even is it an array with `nothing` elements.

### Example Calls
    
Generate an example 2D mesh.
```
x_coords = LinRange(-2, 2, 50)
y_coords = [1 2 4 8 16] # exponentially spaced in y
coords = [x_coords y_coords]
```

Declare two curves, one that does not require parameters and another that does.
```
circle(s) = [cos(2*pi*s), sin(2*pi*s)]                             # Circle of radius r=1 centered at (x,y) = (0,0); this curve does not need parameters
ellipse(s, params) = [param.rx*cos(2*pi*s), param.ry*sin(2*pi*s)]  # Ellispe centered at (0,0) with user-provided radii; this curve takes parameters
```

Find intersections for a single curve:
```
ds = 0.001 # check every 0.36deg
arc_tol, single_tol = 1e-8, 1e-8

intersections = find_path_intersections_single(coords, circle, ds, arc_tol, single_tol) # do not need to pass any parameters
```

Find intersections for multiple curves:
```
ds = [0.001, 0.0005, 0.01] # Can use different step sizes for each curve
arc_tol = [1e-8, 1e-8, 1e-8] # use the same stopping tolerance though
single_tol = [1e-8, 1e-8, 1e-8]
curve_params = [nothing, (0.5, 1), (0.6, 0.2)]
intersections_by_curve = find_path_intersections(coords [circle, ellipse, ellipse], ds, arc_tol, single_tol, curve_params)
```

### Tolerance Defitions
- `arc_tol`
The main tolerance used to find intersections. This tolerance is applied to the Euclidean distance between the upper and lower bounds of the secant-bisection method at each iteration. The method stops when the Euclidean distance between these points is less than or equal to `arc_tol`.

<img src="/figures/tolerance_arcTol.png" width="500">


- `single_tol`
The secondary tolerance for corner intersections. This tolerance is used to determine if another dimension participated in previously found intersection (i.e. to determine if the intersection went through a corner). It is applied to the distance between the previously estimated intersection point and the nearest mesh edge in the dimension being checked. It is, with the exception of folded intersections, a weaker condition than `arc_tol`.

<img src="/figures/tolerances_corners.png" width="900">

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
