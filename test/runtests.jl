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
    circleParam1 = (; r=0.5, x0=0, y0=0.2)
    intersections = find_mesh_intersections_single(coords, circle, ds, arc_tol, single_tol, circleParam1 )

    @test length(intersections) == 4
    @test intersections.s[1] == 0.25
    @test intersections.s[3] == 0.75

    # 2. Simple, single-dim intersections with start/endpoint on a boundary
    # Note: start and end point are both counted, so their common intersection
    # gets counted twice
    circleParam2 = (; r=0.5, x0=0, y0=0)
    intersections = find_mesh_intersections_single(coords, circle, ds, arc_tol, single_tol, circleParam2 )
    @test length(intersections) == 5
    @test intersections.s[1] == 0
    @test intersections.s[2] == 0.25
    @test intersections.s[3] == 0.5
    @test intersections.s[4] == 0.75
    @test intersections.s[5] == 1.0

    # 3. Simple, multiple dim intersections
    epsilon = 0.001
    x_coords = [-1, 0, 0.5*sqrt(2)/2-epsilon, 1]
    y_coords = [-1, 0, 0.5*sqrt(2)/2-epsilon, 1]
    ds = 1/100 # = 1deg
    coords = [x_coords, y_coords]
    circleParams3 = (; r=0.5, x0=0, y0=0)
    
    intersections = find_mesh_intersections_single(coords, circle, ds, arc_tol, single_tol, circleParam2 )
    @test length(intersections) == 9
    @test intersections.s[1] == 0
    @test intersections.s[4] == 0.25
    @test intersections.s[6] == 0.5
    @test intersections.s[7] == 0.75
    @test intersections.s[9] == 1.0


    # 4. Corners, x first, then y

    # 5. Corners, y first, then x 

    # 6. Curve is completely outside the domain

    # 7. Curve that starts outside the domain then enters

    # 8. Curve that starts inside the domain then exits
end
