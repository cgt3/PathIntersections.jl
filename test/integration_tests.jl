using SparseArrays

@testset "Integration Tests" begin
    @testset "Test quadrature on a circle" begin
        circle = PresetGeometries.Circle(R=0.5, x0=0, y0=0)
        all_curves = [circle]

        x_coords = LinRange(-1,1, 6)
        y_coords = LinRange(-1,1, 6)
        coords = [x_coords, y_coords]

        ref_quad = legendre(5)
        regions, cutcell_indices, cutcell_quad, cutcells = get_cutcell_nodes(coords, all_curves, ref_quad)
    
        # Check that the excluded elements are correct
        excluded_elements = findall(regions .< 0)
        @test length(excluded_elements) == 1
        @test excluded_elements[1] == CartesianIndex(3,3) 
        
        # Check that the cut cells are correct
        cut_elements = findall(regions .> 0)
        @test length(cut_elements) == 8
        @test nnz(cutcell_indices) == 8
        @test length(cutcells) == 8

        # Check that the quadrature is correct
        quad_pts = cutcell_quad[:pts]
        quad_wts = cutcell_quad[:wts]
        quad_n   = cutcell_quad[:n]

        sum_perimeter = 0
        for e = 1:length(quad_wts)
            for f = 1:length(quad_wts[e])
                if typeof(cutcells[e].subcurves[f]) <: PresetGeometries.Circle
                    sum_perimeter += sum(quad_wts[e][f])
                end
            end
        end
        @test abs(sum_perimeter - pi) < TESTING_TOL
    end
end