using PathIntersections
using Test

using Revise

@testset "PathIntersections.jl: find_mesh_intersections_single" begin
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
    circleParam = (; r=0.5, x0=0, y0=0)
    intersections = find_mesh_intersections_single(coords, circle, ds, arc_tol, single_tol, circleParam )
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
    coords = [x_coords, y_coords]
    circleParams = (; r=0.5, x0=0, y0=0)
    
    intersections = find_mesh_intersections_single(coords, circle, ds, arc_tol, single_tol, circleParam )
    @test length(intersections) == 9
    @test intersections.s[1] == 0
    @test intersections.s[4] == 0.25
    @test intersections.s[6] == 0.5
    @test intersections.s[7] == 0.75
    @test intersections.s[9] == 1.0


    # 4. Corners, x first, then y
    x_coords = [-1, 0, 0.5*sqrt(2)/2, 1]
    y_coords = [-1, 0, 0.5*sqrt(2)/2, 1]
    coords = [x_coords, y_coords]
    circleParams = (; r=0.5, x0=0, y0=0)
    
    intersections = find_mesh_intersections_single(coords, circle, ds, arc_tol, single_tol, circleParam )
    @test length(intersections) == 8
    @test intersections.s[1] == 0
    @test intersections.s[2] == 0.125
    @test intersections.s[3] == 0.25
    @test intersections.s[4] == 0.375
    @test intersections.s[5] == 0.5
    @test intersections.s[6] == 0.75
    @test intersections.s[7] == 0.875
    @test intersections.s[8] == 1.0

    # 6. Curve is completely outside the domain
    x_coords = LinRange(-1,1,3)
    y_coords = LinRange(-1,1,3)
    coords = [x_coords, y_coords]
    circleParam = (; r=0.5, x0=2, y0=2)
    
    intersections = find_mesh_intersections_single(coords, circle, ds, arc_tol, single_tol, circleParam )
    @test length(intersections) == 0

    # 7. Curve that starts outside the domain then enters
    circleParam = (; r=0.5, x0=1, y0=0)
    intersections = find_mesh_intersections_single(coords, circle, ds, arc_tol, single_tol, circleParam )
    @test length(intersections) == 3
    @test intersections.s[1] == 0.25
    @test intersections.s[2] == 0.5
    @test intersections.s[3] == 0.75

    # 8. Curve that starts inside the domain then exits
    circleParam = (; r=0.5, x0=-1, y0=0)
    # Note the end point and start point are the same so it gets counted twice
    intersections = find_mesh_intersections_single(coords, circle, ds, arc_tol, single_tol, circleParam )
    @test length(intersections) == 4
    @test intersections.s[1] == 0
    @test intersections.s[2] == 0.25
    @test intersections.s[3] == 0.75
    @test intersections.s[4] == 1.0

    # 9. Check passing no parameters for the curve
    circle_noParams(s) = [0.5*cos(2*pi*s), 0.5*sin(2*pi*s)]
    intersections = find_mesh_intersections_single(coords, circle_noParams, ds, arc_tol, single_tol)
    @test length(intersections) == 5
    @test intersections.s[1] == 0
    @test intersections.s[2] == 0.25
    @test intersections.s[3] == 0.5
    @test intersections.s[4] == 0.75
    @test intersections.s[5] == 1.0

    # 10. Check a circle with a negative orientation
    circle_noParams(s) = [0.5*cos(-2*pi*s), 0.5*sin(-2*pi*s)]
    intersections = find_mesh_intersections_single(coords, circle_noParams, ds, arc_tol, single_tol)
    @test length(intersections) == 5
    @test intersections.s[1] == 0
    @test intersections.s[2] == 0.25
    @test intersections.s[3] == 0.5
    @test intersections.s[4] == 0.75
    @test intersections.s[5] == 1.0

end

@testset "PathIntersections.jl: find_mesh_intersections" begin
    x_coords = LinRange(-1,1,3)
    y_coords = LinRange(-1,1,3)
    coords = [x_coords, y_coords]

    # 1. Check processing multiple curves at once
    circleParams = (; r=0.5, x0=1, y0=0)
    circle(s, params) = [params.r*cos(2*pi*s) + params.x0, params.r*sin(2*pi*s) + params.y0]
    circle_noParams(s) = [0.5*cos(2*pi*s), 0.5*sin(2*pi*s)]
    
    curves = [circle, circle_noParams]
    ds_array = [1/100, 1/100]
    arc_tol_array = [1e-8, 1e-8]
    single_tol_array = [1e-8, 1e-8]
    curveParams = [circleParams, nothing]

    intersections_by_curve = find_mesh_intersections(coords, curves, ds_array, arc_tol_array, single_tol_array, curveParams)
    @test length(intersections_by_curve) == 2
    @test length(intersections_by_curve[1]) == 3
    @test length(intersections_by_curve[2]) == 5
end
