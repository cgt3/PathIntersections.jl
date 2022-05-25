@testset "PiecewiseFunctions" begin
    @testset "Invalid arguments" begin
        # number of stop_pts doesn't match # of subfunctions
        stop_pts = (0, 0.5, 1)
        subcurves = ( x->x, x->1, x-> 3*x )
        sub_bounds = ( (0,1), (1,2) )
        @test_throws ErrorException PiecewiseFunction(stop_pts, subcurves, sub_bounds)
        
        # number of stop_pts doesn't match # of bounds
        subcurves = ( x->x, x->1)
        sub_bounds = ( (0,1), )
        @test_throws ErrorException PiecewiseFunction(stop_pts, subcurves, sub_bounds)
    end

    # Discontinuous
    @testset "Discontinuous function" begin
        subfunctions = (x->x, x->-x)
        stop_pts = (-1, 0, 1)
        sub_bounds = ((-2, 0), (4, 5))
        func = PiecewiseFunction(stop_pts, subfunctions, sub_bounds)
        @test func.is_continuous == false
        @test func(-1) == -2
        @test func(0) == -4
        @test func(1) == -5
    end

    # Continuous
    @testset "Continuous function" begin
        subfunctions = (x->x, x->-x)
        stop_pts = (-1, 0, 1)
        sub_bounds = ((-2, 0), (0, 2))
        func = PiecewiseFunction(stop_pts, subfunctions, sub_bounds)
        @test func.is_continuous == true
        @test func(-1) == -2
        @test func(0) == 0
        @test func(1) == -2
    end

    
    # Not defined at a stop point
    @testset "Point discontinuities" begin
        subfunctions = (x->x, x-> (x == 0) ? NaN : -x)
        stop_pts = (-1, 0, 1)
        sub_bounds = ((-1, 0), (0, 1))
        func = PiecewiseFunction(stop_pts, subfunctions, sub_bounds)
        @test func.is_continuous == false
        @test func(-1) == -1
        @test isnan(func(0))
        @test func(1) == -1

        
        subfunctions = (x->x, x-> (x == 1) ? NaN : -x)
        stop_pts = (-1, 0, 1)
        sub_bounds = ((-1, 0), (0, 1))
        func = PiecewiseFunction(stop_pts, subfunctions, sub_bounds)
        @test func.is_continuous == false
        @test func(-1) == -1
        @test func(0) == 0
        @test isnan(func(1))
    end

    
    # Default sub-bounds
    @testset "Default sub-bounds" begin
        subfunctions = (x->x, x->-x)
        stop_pts = (-1, 0, 1)
        func = PiecewiseFunction(stop_pts, subfunctions)
        @test func.is_continuous == true
        @test func.sub_bounds == [(-1,0), (0,1)]
        @test func(-1) == -1
        @test func(0) == 0
        @test func(1) == -1
    end
end # testset: PiecewiseFunctions


@testset "PiecewiseCurves" begin
    @testset "Invalid arguments" begin
        # Bad start point
        stop_pts = (0.5, 1)
        subcurves = ( s->(2,s), )
        sub_bounds = ((0,0.5),)
        @test_throws ErrorException PiecewiseCurve(stop_pts, subcurves, sub_bounds)

        # Bad end point
        stop_pts = (0, 0.5)
        @test_throws ErrorException PiecewiseCurve(stop_pts, subcurves, sub_bounds)

        # Number of stop points doesn't match number of curves
        stop_pts = (0, 0.5, 1)
        subcurves = (s->(2,s), s->(3,s))
        sub_bounds = ((0,0.5),)
        @test_throws ErrorException PiecewiseCurve(stop_pts, subcurves, sub_bounds)

        # Number of stop points doesn't match number of bounds
        stop_pts = (0, 0.5, 1) 
        sub_bounds = ((0,0.5), (0.5, 1), (1, 1.5))
        @test_throws ErrorException PiecewiseCurve(stop_pts, subcurves, sub_bounds)
    end

    # Discontinuous and not closed
    @testset "Discontinuous and not closed curve" begin
        stop_pts = (0, 0.5, 1)
        subcurves = (s->(2,s), s->(3,s))
        sub_bounds = ((0,0.5), (0.5, 1))
        curve = PiecewiseCurve(stop_pts, subcurves, sub_bounds)
        @test curve.is_continuous == false
        @test curve.is_closed == false
        @test curve(0) == (2,0)
        @test curve(0.5) == (3, 0.5)
        @test curve(1) == (3, 1)
    end
    
    # Continuous and closed
    @testset "Continuous and closed curve" begin
        stop_pts = (0, 0.5, 1)
        subcurves = (s->(cos(2*pi*s), sin(2*pi*s)), s->(s, 0))
        sub_bounds = ((0,0.5), (-1, 1))
        curve = PiecewiseCurve(stop_pts, subcurves, sub_bounds)
        @test curve.is_continuous == true
        @test curve.is_closed == true
        @test curve(0) == (1,0)
        @test curve(0.5) == (-1, 0)
        @test curve(0.75) == (0,0)
        @test curve(1) == (1,0)
    end
    
    # Curves wrap inputs to s in [0,1]
    @testset "Curves wraps s to [0,1]" begin
        stop_pts = (0, 1)
        subcurves = (s->s,)
        curve = PiecewiseCurve(stop_pts, subcurves)
        @test curve(1) == curve(2)
        @test curve(-1) == curve(0)
    end

    # ds - arc length based
    @testset "ds: arc length based" begin
        stop_pts = [0, 0.5, 1]
        # has total length 2, but the second segment moves twice as fast
        subcurves = [s->[1,s], s->[2,2*s]]
        sub_bounds = [[0,1], [1, 1.5]]
        curve = PiecewiseCurve(stop_pts, subcurves, sub_bounds)
        @test curve.is_continuous == false
        @test curve.is_closed == false

        ds_func = ds_by_arc_length(curve, 0.1)
        @test abs(ds_func(0) - 0.1) < TESTING_TOL
        @test abs(ds_func(1) - 0.05) < TESTING_TOL
    end

    # ds - num steps
    @testset "ds: num step based" begin
        stop_pts = [0, 0.5, 1]
        # has total length 3, with the second segment having length 2
        subcurves = [s->[1,s], s->[2,2*s]]
        sub_bounds = [[0,1], [1, 2]]
        curve = PiecewiseCurve(stop_pts, subcurves, sub_bounds)
        @test curve.is_continuous == false
        @test curve.is_closed == false

        ds_func = ds_by_num_steps(curve, 100)
        @test abs(ds_func(0) - 0.5/33) < TESTING_TOL
        @test abs(ds_func(1) - 0.5/67) < TESTING_TOL
    end
end