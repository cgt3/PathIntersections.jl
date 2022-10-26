using PathIntersections
using Test

using GaussQuadrature
using LinearAlgebra

TESTING_TOL = 1e-12

include("test_find_mesh_intersections.jl")
include("test_PiecewiseCurves.jl")
include("test_PresetGeometries.jl")
include("test_map_line_quadrature.jl")   

include("integration_tests.jl")
include("is_contained_tests.jl")