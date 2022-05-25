using PathIntersections
using Test

using LinearAlgebra
using Revise

TESTING_TOL = 1e-12

@testset "PathIntersections.jl" begin
    include("test_find_mesh_intersections.jl")
    include("test_PiecewiseCurves.jl")
    include("test_PresetGeometries.jl")
    include("test_map_line_quadrature.jl")    
end