@testset "PresetGeometries" begin
    # Unit circle
    @testset "Default circle" begin
        circle = PresetGeometries.Circle()
        @test norm(circle(0)    .- (1, 0)) < TESTING_TOL
        @test norm(circle(0.25) .- (0, 1)) < TESTING_TOL
        @test norm(circle(0.5)  .- (-1, 0)) < TESTING_TOL
        @test norm(circle(0.75) .- (0, -1)) < TESTING_TOL
        @test norm(circle(1.0)  .- (1, 0)) < TESTING_TOL
    end

    # Default pacman
    @testset "Default pacman" begin
        circle = PresetGeometries.Circle()
        pacman = PresetGeometries.Pacman()
        @test length(pacman.func.stop_pts) == 4
        @test norm(pacman(pacman.func.stop_pts[2]) .- circle(0.125)) < TESTING_TOL
        @test norm(pacman(pacman.func.stop_pts[3]) .- circle(0.875)) < TESTING_TOL
        @test pacman.func.is_continuous == true
        @test pacman.func.is_closed == true
        @test norm(pacman(0)   .- (0, 0)) < TESTING_TOL
        @test norm(pacman(0.5) .- (-1, 0)) < TESTING_TOL
        @test norm(pacman(1)   .- (0, 0)) < TESTING_TOL
        @test abs(pacman(0.1)[1]  .- pacman(0.1)[2]) < TESTING_TOL
        @test abs(pacman(-0.1)[1] .+ pacman(-0.1)[2]) < TESTING_TOL
    end

    
    # Pacman with non-standard angles
    @testset "Pacman: non-standard angles" begin
        # Negative angles
        pacman = PresetGeometries.Pacman(first_jaw=-7*pi/4, second_jaw=-pi/4) # same as 45->(-45)
        @test length(pacman.func.stop_pts) == 4
        @test abs(pacman.first_jaw - pi/4) < TESTING_TOL
        @test abs(pacman.second_jaw - 7*pi/4) < TESTING_TOL

        # Angles over 2pi
        pacman = PresetGeometries.Pacman(first_jaw=9*pi/4, second_jaw=15pi/4) # same as 45->(-45)
        @test length(pacman.func.stop_pts) == 4
        @test abs(pacman.first_jaw - pi/4) < TESTING_TOL
        @test abs(pacman.second_jaw - 7*pi/4) < TESTING_TOL
    end

    
    # "Pacman" can also be a pizza slice
    @testset "Pizza slice Pacman" begin
        # Pizza slice via negative orientation
        pacman = PresetGeometries.Pacman(first_jaw=pi/4, second_jaw=-pi/4, orientation=-1)
        @test length(pacman.func.stop_pts) == 4
        @test abs(pacman.first_jaw - pi/4) < TESTING_TOL
        @test abs(pacman.second_jaw - 7*pi/4) < TESTING_TOL
        @test norm(pacman(0.5) .- (1,0)) < TESTING_TOL

        # Pizza slice via positive orientation
        pacman = PresetGeometries.Pacman(first_jaw=pi/2, second_jaw=0) 
        @test length(pacman.func.stop_pts) == 4
        @test abs(pacman.first_jaw - pi/2) < TESTING_TOL
        @test abs(pacman.second_jaw - 0) < TESTING_TOL
        @test norm(pacman(0.5) .+ (sqrt(2)/2, sqrt(2)/2)) < TESTING_TOL
    end

    # Ellipse
    @testset "Ellipse, default" begin
        ellipse = PresetGeometries.Ellipse(Rx=5, Ry=3, x0=0, y0=5)
        @test norm(ellipse(0)    .- (5, 5)) < TESTING_TOL
        @test norm(ellipse(0.25) .- (0, 8)) < TESTING_TOL
        @test norm(ellipse(0.5)  .- (-5, 5)) < TESTING_TOL
        @test norm(ellipse(0.75) .- (0, 2)) < TESTING_TOL
        @test norm(ellipse(1)    .- (5, 5)) < TESTING_TOL
    end

    @testset "Ellipse, rotated" begin
        ellipse = PresetGeometries.Ellipse(Rx=5, Ry=3, theta0=pi/2)
        @test norm(ellipse(0)    .- (0,5)) < TESTING_TOL
        @test norm(ellipse(0.25) .- (-3,0)) < TESTING_TOL
        @test norm(ellipse(0.5)  .- (0,-5)) < TESTING_TOL
        @test norm(ellipse(0.75) .- (3,0)) < TESTING_TOL
        @test norm(ellipse(1)    .- (0,5)) < TESTING_TOL
    end
end # testset: PresetGeometries