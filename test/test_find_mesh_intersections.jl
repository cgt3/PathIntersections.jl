
@testset "find_mesh_intersections" begin
    x_coords = LinRange(-1,1,3)
    y_coords = LinRange(-1,1,3)
    coords = [x_coords, y_coords]
    
    ds = 1/100 # check every 3.6deg
    arc_tol = 1e-8
    corner_tol = 1e-8

    circle(s, param) = SVector(param.r*cos(2*pi*s) + param.x0, param.r*sin(2*pi*s) + param.y0)
    circle_noParams(s) = SVector(0.5*cos(2*pi*s), 0.5*sin(2*pi*s))

    @testset "Simple, single-dim" begin
        circle1(s) = circle(s, (; r=0.5, x0=0, y0=0.2))
        intersections = find_mesh_intersections(coords, circle1, ds, arc_tol, corner_tol)

        @test length(intersections) == 4
        @test intersections[1].s == 0.25
        @test intersections[3].s == 0.75
    end

    @testset "Simple, single-dim, closed-list enforced" begin
        circle1(s) = circle(s, (; r=0.5, x0=0, y0=0.2))
        intersections = find_mesh_intersections(coords, circle1, ds, arc_tol, corner_tol, closed_list=true)

        @test length(intersections) == 5
        @test intersections[1].s == 0.25
        @test intersections[3].s == 0.75
        @test intersections[5].s == 0.25
    end

    
    @testset "Piecewise curve with stop points" begin
        rect = PresetGeometries.Rectangle(Lx=0.5, Ly=0.5)
        intersections = find_mesh_intersections(coords, rect, ds, arc_tol, corner_tol, closed_list=true)

        @test length(intersections) == 9
        # The corners
        @test intersections[1].s == 0
        @test intersections[3].s == 0.25
        @test intersections[5].s == 0.5
        @test intersections[7].s == 0.75
        @test intersections[9].s == 1

        # The intersections
        @test intersections[2].s == 0.125
        @test intersections[4].s == 0.375
        @test intersections[6].s == 0.625
        @test intersections[8].s == 0.875
    end

    # Note: start and end point are both counted, so their common intersection
    # gets counted twice
    @testset "Simple, single-dim; start/end pt on boundary" begin
        circle2(s) = circle(s, (; r=0.5, x0=0, y0=0))
        intersections = find_mesh_intersections(coords, circle2, ds, arc_tol, corner_tol)
        @test length(intersections) == 5
        @test intersections[1].s == 0
        @test intersections[2].s == 0.25
        @test intersections[3].s == 0.5
        @test intersections[4].s == 0.75
        @test intersections[5].s == 1.0
    end

    # 3. Simple, multiple dim intersections
    @testset "Simple, multiple-dim" begin
        epsilon = 0.001
        x_coords = [-1, 0, 0.5*sqrt(2)/2-epsilon, 1]
        y_coords = [-1, 0, 0.5*sqrt(2)/2-epsilon, 1]
        coords = [x_coords, y_coords]
        circle3(s) = circle(s, (; r=0.5, x0=0, y0=0))
        
        intersections = find_mesh_intersections(coords, circle3, ds, arc_tol, corner_tol)
        @test length(intersections) == 9
        @test intersections[1].s == 0
        @test intersections[4].s == 0.25
        @test intersections[6].s == 0.5
        @test intersections[7].s == 0.75
        @test intersections[9].s == 1.0
    end

    # 4. Corners, x first, then y
    @testset "Corners" begin
        x_coords = [-1, 0, 0.5*sqrt(2)/2, 1]
        y_coords = [-1, 0, 0.5*sqrt(2)/2, 1]
        coords = [x_coords, y_coords]
        circle4(s) = circle(s, (; r=0.5, x0=0, y0=0))
        
        intersections = find_mesh_intersections(coords, circle4, ds, arc_tol, corner_tol)
        @test length(intersections) == 8
        @test intersections[1].s == 0
        @test intersections[2].s == 0.125
        @test intersections[3].s == 0.25
        @test intersections[4].s == 0.375
        @test intersections[5].s == 0.5
        @test intersections[6].s == 0.75
        @test intersections[7].s == 0.875
        @test intersections[8].s == 1.0
    end

    # 5. Curve is completely outside the domain
    @testset "Curve out of bounds" begin
        x_coords = LinRange(-1,1,3)
        y_coords = LinRange(-1,1,3)
        coords = [x_coords, y_coords]
        circle5(s) = circle(s, (; r=0.5, x0=2, y0=2))
        
        intersections = find_mesh_intersections(coords, circle5, ds, arc_tol, corner_tol)
        @test length(intersections) == 0
    end

    # 6. Curve that starts outside the domain then enters
    @testset "Curve starts outside the domain" begin
        circle6(s) = circle(s, (; r=0.5, x0=1, y0=0))
        intersections = find_mesh_intersections(coords, circle6, ds, arc_tol, corner_tol)
        @test length(intersections) == 3
        @test intersections[1].s == 0.25
        @test intersections[2].s == 0.5
        @test intersections[3].s == 0.75
    end

    # 7. Curve that starts inside the domain then exits
    @testset "Curve start inside domain then leaves" begin
        circle7(s) = circle(s, (; r=0.5, x0=-1, y0=0))
        # Note the end point and start point are the same so it gets counted twice
        intersections = find_mesh_intersections(coords, circle7, ds, arc_tol, corner_tol)
        @test length(intersections) == 4
        @test intersections[1].s == 0
        @test intersections[2].s == 0.25
        @test intersections[3].s == 0.75
        @test intersections[4].s == 1.0
    end

    # 9. Complete code coverage by tripping intersections on lower boundaries
    @testset "Lower bound intersections" begin            
        circle_neg(s) = SVector(0.5*cos(-2*pi*s), 0.5*sin(-2*pi*s))
        
        x_coords = [-1 0 1]
        y_coords = [-1 -0.5 circle_neg(2*ds)[2] 0 1]
        coords = [x_coords, y_coords]
        intersections = find_mesh_intersections(coords, circle_neg, ds, arc_tol, corner_tol)
        @test length(intersections) == 7
        @test intersections[1].s == 0
        @test intersections[2].s == 2*ds
        @test intersections[3].s == 0.25
        @test intersections[5].s == 0.5
        @test intersections[6].s == 0.75
        @test intersections[7].s == 1.0
    end

    
    # 10. Provide only the step size and allow the tolerances to default
    @testset "Default tolerances" begin            
        circle10(s) = circle(s, (; r=0.5, x0=0, y0=0))
        x_coords = [-1.0 0.0 1.0]
        y_coords = [-1.0 0.0 1.0]
        coords = [x_coords, y_coords]

        intersections = find_mesh_intersections(coords, circle10, ds)
        @test length(intersections) == 5
        @test intersections[1].s == 0
        @test intersections[2].s == 0.25
        @test intersections[3].s == 0.5
        @test intersections[4].s == 0.75
        @test intersections[5].s == 1.0

    end

    # 11. Use default step size
    @testset "Default step size and tolerances" begin            
        circle11(s) = circle(s, (; r=0.5, x0=0, y0=0))
        x_coords = [-1.0 0.0 1.0]
        y_coords = [-1.0 0.0 1.0]
        coords = [x_coords, y_coords]

        intersections = find_mesh_intersections(coords, circle11)
        @test length(intersections) == 5
        @test intersections[1].s == 0
        @test intersections[2].s == 0.25
        @test intersections[3].s == 0.5
        @test intersections[4].s == 0.75
        @test intersections[5].s == 1.0
    end

    
    # 12. Using a function for the step size
    @testset "Providing ds as a function; single curve" begin
        # A curve where the arc length per constant step in s varies        
        ellipse(s) = SVector(0.5*cos(2*pi*s), 0.05*sin(2*pi*s))
        ds_func(s) = ds + ds*cos(2*pi*s)^2

        x_coords = [-1.0 0.0 1.0]
        y_coords = [-1.0 0.0 1.0]
        coords = [x_coords, y_coords]

        intersections = find_mesh_intersections(coords, ellipse, ds_func, arc_tol, corner_tol)
        @test length(intersections) == 5
        @test intersections[1].s == 0
        @test intersections[2].s == 0.25
        @test intersections[3].s == 0.5
        @test intersections[4].s == 0.75
        @test intersections[5].s == 1.0

    end
    
    # Tests sending multiple curves at once ------------------------------------------------
    @testset "Array of curves" begin
        x_coords = LinRange(-1,1,3)
        y_coords = LinRange(-1,1,3)
        coords = [x_coords, y_coords]

        # 1. Check processing multiple curves at once
        circle1(s) = circle(s, (; r=0.5, x0=1, y0=0))

        curves = [circle1, circle_noParams]
        ds_array = [1/100, 1/100]
        arc_tol_array = [1e-8, 1e-8]
        corner_tol_array = [1e-8, 1e-8]

        intersections_by_curve = find_mesh_intersections(coords, curves, ds_array,
            arc_tol_array, corner_tol_array)
        @test length(intersections_by_curve) == 2
        @test length(intersections_by_curve[1]) == 3
        @test length(intersections_by_curve[2]) == 5
    end

    @testset "Tuple of curves" begin
        x_coords = LinRange(-1,1,3)
        y_coords = LinRange(-1,1,3)
        coords = [x_coords, y_coords]

        # 1. Check processing multiple curves at once
        circle1(s) = circle(s, (; r=0.5, x0=1, y0=0))

        curves = (circle1, circle_noParams)
        ds_array = (1/100, 1/100)
        arc_tol_array = (1e-8, 1e-8)
        corner_tol_array = (1e-8, 1e-8)

        intersections_by_curve = find_mesh_intersections(coords, curves, ds_array,
            arc_tol_array, corner_tol_array)
        @test length(intersections_by_curve) == 2
        @test length(intersections_by_curve[1]) == 3
        @test length(intersections_by_curve[2]) == 5
    end
    
    @testset "Multiple curves, scalar step size and tolerances" begin
        x_coords = LinRange(-1,1,3)
        y_coords = LinRange(-1,1,3)
        coords = [x_coords, y_coords]

        # 1. Check processing multiple curves at once
        circle1(s) = circle(s, (; r=0.5, x0=1, y0=0))

        curves = (circle1, circle_noParams)
        ds = 1/100
        arc_tol = 1e-8
        corner_tol = 1e-8

        intersections_by_curve = find_mesh_intersections(coords, curves, ds, arc_tol, corner_tol)
        @test length(intersections_by_curve) == 2
        @test length(intersections_by_curve[1]) == 3
        @test length(intersections_by_curve[2]) == 5
    end
    
    @testset "Multiple curves, default step size and tolerances" begin
        x_coords = LinRange(-1,1,3)
        y_coords = LinRange(-1,1,3)
        coords = [x_coords, y_coords]

        # 1. Check processing multiple curves at once
        circle1(s) = circle(s, (; r=0.5, x0=1, y0=0))

        curves = (circle1, circle_noParams)
        intersections_by_curve = find_mesh_intersections(coords, curves)
        @test length(intersections_by_curve) == 2
        @test length(intersections_by_curve[1]) == 3
        @test length(intersections_by_curve[2]) == 5
    end

    
    @testset "Providing ds as a function; multiple curves" begin
        x_coords = LinRange(-1,1,3)
        y_coords = LinRange(-1,1,3)
        coords = [x_coords, y_coords]

        # 1. Check processing multiple curves at once
        circle1(s) = circle(s, (; r=0.5, x0=1, y0=0))
        ds_func(s) = ds + ds*cos(2*pi*s)^2

        curves = (circle1, circle_noParams)
        ds = 1/100
        arc_tol = 1e-8
        corner_tol = 1e-8

        intersections_by_curve = find_mesh_intersections(coords, curves, ds_func, arc_tol, corner_tol)
        @test length(intersections_by_curve) == 2
        @test length(intersections_by_curve[1]) == 3
        @test length(intersections_by_curve[2]) == 5
    end
end # testset: find_mesh_intersections