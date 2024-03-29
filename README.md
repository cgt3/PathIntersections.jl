# PathIntersections

[![Build status](https://github.com/cgt3/PathIntersections.jl/workflows/CI/badge.svg)](https://github.com/cgt3/PathIntersections.jl/actions)
[![Build Status](https://ci.appveyor.com/api/projects/status/github/cgt3/PathIntersections.jl?svg=true)](https://ci.appveyor.com/project/cgt3/PathIntersections-jl)
[![Coverage](https://codecov.io/gh/cgt3/PathIntersections.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/cgt3/PathIntersections.jl)

## Introduction
This package locates intersections between parameterized curves (paths) and Cartesian meshes in arbitrary dimension. The package locates intersections by walking along the curve using the user provided step size, `ds`, and sensing for boundary crossings. When a boundary is crossed, a hybrid bracketed secant-bisection method is used to determine the location of the intersection within the user-specified tolerances `arc_tol` and `corner_tol`. 

While the use of a bracketed method is robust in many situations, it does not allow all non-simple intersections to be found; cases that are currently unsupported or may result in unsatisfactory accuracy are noted below.

## Usage and Examples
The package provides the function `find_mesh_intersections` for finding intersections between a mesh and a single curve or a set of curves.

### Return values:
- **Single curve**: An array of type `MeshIntersection`'s where the struct `MeshIntersections` is defined:
    ```
    mutable struct MeshIntersection
        s::Real                      # The s-value where the intersection occurred
        pt::AbstractArray{Real}      # The approximate intersection point
        dim::Vector{Bool}            # Whether each dimension was involved in the intersection
        I::AbstractArray{Real}       # The indices of the closest lower mesh bound or the       
                                     # boundaries involved in the intersection
    end
    ```

- **Multiple curves**: An array of arrays of type `MeshIntersection`'s. There is one array per curve with indexing: `intersections_by_curve[curve index][intersection index]`.

### Arguments:
- `coords`: A Cartesian mesh, which is to be provided in the form of an array of arrays with indexing `coords[dimension][edge index]`. Note the mesh does not need to be uniform.

- `curve`/`curves` : The curve/curves to be check against the mesh. The package expects curves to be callable using `pt = curve(s)` where `pt` is indexed `pt[dimension]`. For an array/tuple of curves the expected indexing/signature is `pt = curves[curve index](s)`.

- `ds` : optional argument, the/an array of step size(s) use for walking along the curve. Arrays are to be indexed `ds[curve index]`. Defaults to the minimum of 1/1000 (note `s` is limited to [0,1]) or a value such that `norm(curve(s+ds) - curve(s))` is less than the minimum mesh spacing.

- `arc_tol` : optional argument, the/an array of arc length tolerance(s) used as the stopping criteria for finding intersections. Arrays are to be indexed `arc_tol[curve index]`. Default is `100*eps(eltype(coords))`.

- `corner_tol`: optional argument, the/an array of single-dimension tolerance(s) used for identifying if other dimensions participated in an intersecton (i.e. whether the intersection occurred at a corner). Arrays are to be indexed `corner_tol[curve index]`. Default is `100*eps(eltype(coords))`.


### Assumptions
1. Each array in the array of mesh coordinates, `coords`, must:
    1. be sorted in increasing order (i.e. `coords[dim][i] < coords[dim][i+1]`),
    2. not contain duplicates, and
    3. contain the start and end points of the domain in that dimension.
2. All curves start at `s = 0` and end at `s = 1` and are defined at every point in between. The end points will be checked regardless of the step size provided.
3. The provided step size is sufficiently small in comparison to the mesh for all intersections to be sensed.
4. Differentiability of the curve: When the step size, `ds`, is defaulted, the derivative of the curve will be evaluated in order to find a suitable step size using automatic differentiation as implemented in `ForwardDiff.jl`. If the curve or its derivative is not defined at certain points the resulting step size may be undesirable. 

### Notes:
1. This package is (currently) only concerned with mesh-curve intersections, not curve-curve intersections.
2. Intersections are returned sorted in increasing order according to their `s` values.
3. Curves may start and arbitrarily enter and leave regions outside the domain defined by the mesh.


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
arc_tol, corner_tol = 1e-8, 1e-8

intersections = find_mesh_intersections(coords, circle, ds, arc_tol, corner_tol)
```

Find intersections for multiple curves:
```
ds = [0.001, 0.0005, 0.01] # Can use different step sizes for each curve
arc_tol = [1e-8, 1e-8, 1e-8] # Can also use all the same
corner_tol = [1e-8, 1e-8, 1e-8]
ellipse1(s) = ellipse(s, (0.5,1))
ellipse2(s) = ellipse(s, (0.6, 0.2))
intersections_by_curve = find_mesh_intersections(coords, [circle, ellipse1, ellipse2], ds, arc_tol, corner_tol)
```

### Tolerance Definitions
- `arc_tol`
The main tolerance used to find intersections. This tolerance is applied to the Euclidean distance between the upper and lower bounds of the secant-bisection method at each iteration. The method stops when the Euclidean distance between these points is less than or equal to `arc_tol`.

<p align="center">
  <img src="/figures/tolerance_arcTol.png" width="800">
</p>



- `corner_tol`
The secondary tolerance for corner intersections. This tolerance is used to determine if another dimension participated in previously found intersection (i.e. to determine if the intersection went through a corner). It is applied to the distance between the previously estimated intersection point and the nearest mesh edge in the dimension being checked. It is, with the exception of folded intersections, a weaker condition than `arc_tol`.

<p align="center">
    <img src="/figures/tolerances_corners.png" width="900">
</p>

## Supported Cases
### Simple, Single-Dimension Intersections

<p align="center">
    <img src="/figures/intersection_simpleSingleDim.png" width="700">
</p>

### Simple, Multiple-Dimension Intersections

<p align="center">
    <img src="/figures/intersection_simpleMultiDim.png" width="700">
</p>

### Corners

<p align="center">
    <img src="/figures/intersection_corner.png" width="700">
</p>

## Unsupported and Degenerate Cases
### Unsupported: Non-simple Intersections

<p align="center">
    <img src="/figures/intersections_nonSimple.png" width="700">
</p>

### Unsupported: Bad Step Size
Note if the step size is too large wrt the mesh, intersections may be missed. Missed intersections will also occur if both the current and next point are on boundaries.
<p align="center">
    <img src="/figures/intersections_badStepSize.png" width="700">
</p>

### Degenerate Input: Folded Intersections

<p align="center">
    <img src="/figures/tolerance_foldedIntersection.png" width="800">
</p>
