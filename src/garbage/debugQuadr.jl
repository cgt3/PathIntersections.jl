using PathIntersections

using GaussQuadrature
using LinearAlgebra
using Revise

# @testset "map_line_quadrature (multiple pts)" begin
    ref_pts, ref_wts = legendre(5) # Legendre-Gauss quadrature
    true_pts_2D = [ 0.5*[ref_pts[i] + 1, ref_pts[i] + 1] for i in 1:length(ref_pts)]

    line(s) = [s, s] # simple linear function
    s_bounds = [(-1, 0), (0, 1)] # Yields 2 intervals, [-1,0], [0,1]
    
    line_pts, line_wts, normals = map_line_quadrature_nonconsecutive((ref_pts, ref_wts), line, s_bounds, normalization=true)
    # @test length(line_wts) == 2
    # @test length(line_pts) == 2
    # @test length(normals) == 2
    
    # @test length(line_wts[1]) == 5
    # @test length(line_pts[1]) == 5
    # @test length(normals[1]) == 5
    
    # @test length(line_wts[2]) == 5
    # @test length(line_pts[2]) == 5
    # @test length(normals[2]) == 5

    # # Interval [-1, 0]
    # for i = 1:length(ref_pts)
    #     diff_pts = abs.(line_pts[1][i] - (-true_pts_2D[end-i+1]) )
    #     @test maximum(diff_pts) < TESTING_TOL
            
    #     diff_normals = abs.(normals[1][i] - sqrt(2)/2*[1, -1])
    #     @test maximum(diff_normals) < TESTING_TOL
    # end

    # # Interval [0, 1]
    # for i = 1:length(ref_pts)
    #     diff_pts = abs.(line_pts[2][i] - true_pts_2D[i])
    #     @test maximum(diff_pts) < TESTING_TOL
            
    #     diff_normals = abs.(normals[2][i] - sqrt(2)/2*[1, -1])
    #     @test maximum(diff_normals) < TESTING_TOL
    # end