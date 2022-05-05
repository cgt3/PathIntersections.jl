using PathIntersections
using Test

using Revise

@testset "PathIntersections.jl" begin
    # Write your tests here.
    intersection = Intersection(0, [1,1], [1,0], [2,3])
    @test intersection.s == 0
    
    x_coords = LinRange(-1,1,3)
    y_coords = LinRange(-1,1,3)
    coords = [x_coords, y_coords]
    
    ds = 1/100 # check every 3.6deg
    arc_tol = 1e-8
    single_tol = 1e-8
    circle(s, param) = [param.r*cos(2*pi*s) + param.x0, param.r*sin(2*pi*s) + param.y0]

    # 1. Simple, single-dim intersections, start and endpoints not on boundaries
    circleParam1 = (; r=0.5, x0=-0.2, y0=0.2)
    intersections = find_mesh_intersections_single(coords, circle, ds, arc_tol, single_tol, circleParam1 )

    @test length(intersections) == 4

    # 2. Simple, single-dim intersections with start/endpoint on a boundary
    # Note: start and end point are both counted, so their common intersection
    # gets counted twice
    circleParam2 = (; r=0.5, x0=0, y0=0)
    intersections = find_mesh_intersections_single(coords, circle, ds, arc_tol, single_tol, circleParam2 )
    print(intersections)
    @test length(intersections) == 5
end
