
# push!(LOAD_PATH, "../src")
#import Pkg; Pkg.add("PathIntersections")
using PathIntersections
using Test

@testset "PathIntersections.jl" begin
    # Write your tests here.
    intersection = Intersection(0, [1,1], [1,0], [2,3])
    @test intersection.s == 0
end
