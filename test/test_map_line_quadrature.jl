using GaussQuadrature

TESTING_TOL = 1e-15

@testset "Quadrature Mapping" begin
    @testset "map_line_quadrature" begin
        # 1. Identity mapping
        @testset "Identity mapping" begin
            ref_pts, ref_wts = legendre(5) # Legendre-Gauss quadrature
            ref_pts_2D = [ [ref_pts[i], 0] for i in 1:length(ref_pts)]

            line(s) = [s, 0*s] # the x-axis
            s_domain = (-1, 1)

            line_pts, line_wts, normals = map_line_quadrature((ref_pts, ref_wts), line, s_domain)
            @test line_wts == ref_wts
            @test length(line_pts) == length(ref_pts)

            for i = 1:length(ref_pts)
                diff = abs.(line_pts[i] - ref_pts_2D[i])
                @test maximum(diff) < TESTING_TOL
                @test normals[i] == [0, -1]
            end
        end

        # 2. Linear mapping to test normalization
        @testset "Linear mapping w normalization" begin
            ref_pts, ref_wts = legendre(5) # Legendre-Gauss quadrature
            true_pts_2D = [ [ref_pts[i], ref_pts[i]] for i in 1:length(ref_pts)]

            line(s) = [s, s] # simple linear function
            s_domain = (-1, 1)

            line_pts, line_wts, normals = map_line_quadrature((ref_pts, ref_wts), line, s_domain, normalization=true)
            @test line_wts == sqrt(2)*ref_wts
            @test length(line_pts) == length(ref_pts)

            for i = 1:length(ref_pts)
                diff_pts = abs.(line_pts[i] - true_pts_2D[i])
                @test maximum(diff_pts) < TESTING_TOL
                
                diff_normals = abs.(normals[i] - sqrt(2)/2*[1, -1])
                @test maximum(diff_normals) < TESTING_TOL
            end
        end
    end # testset: "map_line_quadrature"

    @testset "map_line_quadrature_consecutive" begin
        ref_pts, ref_wts = legendre(5) # Legendre-Gauss quadrature
        true_pts_2D = [ 0.5*[ref_pts[i] + 1, ref_pts[i] + 1] for i in 1:length(ref_pts)]

        line(s) = [s, s] # simple linear function
        stop_pts = [-1, 0, 1] # Yields 2 intervals, [-1,0], [0,1]

        line_pts, line_wts, normals = map_line_quadrature_consecutive((ref_pts, ref_wts), line, stop_pts, normalization_all=true)
        @test length(line_wts) == 2
        @test length(line_pts) == 2
        @test length(normals) == 2
        
        @test length(line_wts[1]) == 5
        @test length(line_pts[1]) == 5
        @test length(normals[1]) == 5
        
        @test length(line_wts[2]) == 5
        @test length(line_pts[2]) == 5
        @test length(normals[2]) == 5

        # Interval [-1, 0]
        for i = 1:length(ref_pts)
            diff_pts = abs.(line_pts[1][i] - (-true_pts_2D[end-i+1]) )
            @test maximum(diff_pts) < TESTING_TOL
                
            diff_normals = abs.(normals[1][i] - sqrt(2)/2*[1, -1])
            @test maximum(diff_normals) < TESTING_TOL
        end

        # Interval [0, 1]
        for i = 1:length(ref_pts)
            diff_pts = abs.(line_pts[2][i] - true_pts_2D[i])
            @test maximum(diff_pts) < TESTING_TOL
                
            diff_normals = abs.(normals[2][i] - sqrt(2)/2*[1, -1])
            @test maximum(diff_normals) < TESTING_TOL
        end
    end # testset: "map_consecutive_line_quadratures

    @testset "map_line_quadrature_multiple" begin
        ref_pts, ref_wts = legendre(5) # Legendre-Gauss quadrature
        true_pts_2D = [ 0.5*[ref_pts[i] + 1, ref_pts[i] + 1] for i in 1:length(ref_pts)]

        line(s) = [s, s] # simple linear function
        s_bounds = [(-1, 0), (0, 1)] # Yields 2 intervals, [-1,0], [0,1]

        line_pts, line_wts, normals = map_line_quadrature_multiple((ref_pts, ref_wts), line, s_bounds, normalization_all=true)
        @test length(line_wts) == 2
        @test length(line_pts) == 2
        @test length(normals) == 2
        
        @test length(line_wts[1]) == 5
        @test length(line_pts[1]) == 5
        @test length(normals[1]) == 5
        
        @test length(line_wts[2]) == 5
        @test length(line_pts[2]) == 5
        @test length(normals[2]) == 5

        # Interval [-1, 0]
        for i = 1:length(ref_pts)
            diff_pts = abs.(line_pts[1][i] - (-true_pts_2D[end-i+1]) )
            @test maximum(diff_pts) < TESTING_TOL
                
            diff_normals = abs.(normals[1][i] - sqrt(2)/2*[1, -1])
            @test maximum(diff_normals) < TESTING_TOL
        end

        # Interval [0, 1]
        for i = 1:length(ref_pts)
            diff_pts = abs.(line_pts[2][i] - true_pts_2D[i])
            @test maximum(diff_pts) < TESTING_TOL
                
            diff_normals = abs.(normals[2][i] - sqrt(2)/2*[1, -1])
            @test maximum(diff_normals) < TESTING_TOL
        end
    end # testset: "map_multiple_line_quadratures
end
