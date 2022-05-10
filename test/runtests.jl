using PathIntersections
using Test

using Revise

@testset "PathIntersections.jl" begin
    @testset "find_mesh_intersections" begin
        x_coords = LinRange(-1,1,3)
        y_coords = LinRange(-1,1,3)
        coords = [x_coords, y_coords]
        
        ds = 1/100 # check every 3.6deg
        arc_tol = 1e-8
        single_tol = 1e-8

        circle(s, param) = [param.r*cos(2*pi*s) + param.x0, param.r*sin(2*pi*s) + param.y0]
        circle_noParams(s) = [0.5*cos(2*pi*s), 0.5*sin(2*pi*s)]

        @testset "Simple, single-dim" begin
            circleParam1 = (; r=0.5, x0=0, y0=0.2)
            intersections = find_mesh_intersections(coords, circle, ds, arc_tol, single_tol, circleParam1 )

            @test length(intersections) == 4
            @test intersections.s[1] == 0.25
            @test intersections.s[3] == 0.75
        end

        # Note: start and end point are both counted, so their common intersection
        # gets counted twice
        @testset "Simple, single-dim; start/end pt on boundary" begin
            circleParam = (; r=0.5, x0=0, y0=0)
            intersections = find_mesh_intersections(coords, circle, ds, arc_tol, single_tol, circleParam )
            @test length(intersections) == 5
            @test intersections.s[1] == 0
            @test intersections.s[2] == 0.25
            @test intersections.s[3] == 0.5
            @test intersections.s[4] == 0.75
            @test intersections.s[5] == 1.0
        end

        # 3. Simple, multiple dim intersections
        @testset "Simple, multiple-dim" begin
            epsilon = 0.001
            x_coords = [-1, 0, 0.5*sqrt(2)/2-epsilon, 1]
            y_coords = [-1, 0, 0.5*sqrt(2)/2-epsilon, 1]
            coords = [x_coords, y_coords]
            circleParam = (; r=0.5, x0=0, y0=0)
            
            intersections = find_mesh_intersections(coords, circle, ds, arc_tol, single_tol, circleParam )
            @test length(intersections) == 9
            @test intersections.s[1] == 0
            @test intersections.s[4] == 0.25
            @test intersections.s[6] == 0.5
            @test intersections.s[7] == 0.75
            @test intersections.s[9] == 1.0
        end

        # 4. Corners, x first, then y
        @testset "Corners" begin
            x_coords = [-1, 0, 0.5*sqrt(2)/2, 1]
            y_coords = [-1, 0, 0.5*sqrt(2)/2, 1]
            coords = [x_coords, y_coords]
            circleParam = (; r=0.5, x0=0, y0=0)
            
            intersections = find_mesh_intersections(coords, circle, ds, arc_tol, single_tol, circleParam )
            @test length(intersections) == 8
            @test intersections.s[1] == 0
            @test intersections.s[2] == 0.125
            @test intersections.s[3] == 0.25
            @test intersections.s[4] == 0.375
            @test intersections.s[5] == 0.5
            @test intersections.s[6] == 0.75
            @test intersections.s[7] == 0.875
            @test intersections.s[8] == 1.0
        end

        # 6. Curve is completely outside the domain
        @testset "Curve out of bounds" begin
            x_coords = LinRange(-1,1,3)
            y_coords = LinRange(-1,1,3)
            coords = [x_coords, y_coords]
            circleParam = (; r=0.5, x0=2, y0=2)
            
            intersections = find_mesh_intersections(coords, circle, ds, arc_tol, single_tol, circleParam )
            @test length(intersections) == 0
        end

        # 7. Curve that starts outside the domain then enters
        @testset "Curve starts outside the domain" begin
            circleParam = (; r=0.5, x0=1, y0=0)
            intersections = find_mesh_intersections(coords, circle, ds, arc_tol, single_tol, circleParam )
            @test length(intersections) == 3
            @test intersections.s[1] == 0.25
            @test intersections.s[2] == 0.5
            @test intersections.s[3] == 0.75
        end

        # 8. Curve that starts inside the domain then exits
        @testset "Curve start inside domain then leaves" begin
            circleParam = (; r=0.5, x0=-1, y0=0)
            # Note the end point and start point are the same so it gets counted twice
            intersections = find_mesh_intersections(coords, circle, ds, arc_tol, single_tol, circleParam )
            @test length(intersections) == 4
            @test intersections.s[1] == 0
            @test intersections.s[2] == 0.25
            @test intersections.s[3] == 0.75
            @test intersections.s[4] == 1.0
        end

        # 9. Check passing no parameters for the curve
        @testset "No params passed to the curve" begin
            intersections = find_mesh_intersections(coords, circle_noParams, ds, arc_tol, single_tol)
            @test length(intersections) == 5
            @test intersections.s[1] == 0
            @test intersections.s[2] == 0.25
            @test intersections.s[3] == 0.5
            @test intersections.s[4] == 0.75
            @test intersections.s[5] == 1.0
        end

        # 10. Complete code coverage by tripping intersections on lower boundaries
        @testset "Lower bound intersections" begin            
            circle_neg(s) = [0.5*cos(-2*pi*s), 0.5*sin(-2*pi*s)]
            
            x_coords = [-1 0 1]
            y_coords = [-1 -0.5 circle_neg(2*ds)[2] 0 1]
            coords = [x_coords, y_coords]
            intersections = find_mesh_intersections(coords, circle_neg, ds, arc_tol, single_tol)
            @test length(intersections) == 7
            @test intersections.s[1] == 0
            @test intersections.s[2] == 2*ds
            @test intersections.s[3] == 0.25
            @test intersections.s[5] == 0.5
            @test intersections.s[6] == 0.75
            @test intersections.s[7] == 1.0

        end
        
        # Tests sending multiple curves at once ------------------------------------------------
        @testset "Array of curves" begin
            x_coords = LinRange(-1,1,3)
            y_coords = LinRange(-1,1,3)
            coords = [x_coords, y_coords]

            # 1. Check processing multiple curves at once
            circleParams = (; r=0.5, x0=1, y0=0)

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
    end # testset: find_mesh_intersections
end