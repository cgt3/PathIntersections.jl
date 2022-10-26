@testset "`is_contained` tests" begin
    circle = PresetGeometries.Circle()
    @test is_contained(circle, [0., 0.])==true

    # TODO: add an `is_contained` test for a general curve
end