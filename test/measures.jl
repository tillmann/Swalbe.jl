@testset "Measurements" begin
    @testset "Wetted liquid solid" begin
        area_size = [1.0]
        maxheight = [1.0]
        height = zeros(5,5)
        height[2:3, 2:3] .= 1.0
        Swalbe.wetted!(area_size, maxheight, height, 1)
        @test area_size[1] == 4.0
        @test maxheight[1] == 1.0
    
        area_size = ones(5)
        height = zeros(5,5)
        maxheight = zeros(5)
        for t in 1:5
            height[t, t] = t
            Swalbe.wetted!(area_size, maxheight, height, t)
        end
        sol1 = [1.0, 2.0, 3.0, 4.0, 5.0]
        
        for t in 1:5
            @test area_size[t] == sol1[t]
            @test maxheight[t] == sol1[t]
        end
    end

    @testset "Position liquid" begin
        height = zeros(5,5)
        fluid = falses(1,25)
        height[2:3, 2:3] .= 1.0
        dummy = falses(5,5)
        Swalbe.fluid_dry!(fluid, dummy, height, 1)

        @test all(fluid[1, 1:6] .== false)
        @test all(fluid[1, 7:8] .== true)
    
        fluid = falses(5,25)
        height = zeros(5,5)
        maxheight = zeros(5)
        # Desired solution
        sol2 = [1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0; 
                1 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0; 
                1 0 0 0 0 0 1 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0; 
                1 0 0 0 0 0 1 0 0 0 0 0 1 0 0 0 0 0 1 0 0 0 0 0 0; 
                1 0 0 0 0 0 1 0 0 0 0 0 1 0 0 0 0 0 1 0 0 0 0 0 1]
        for t in 1:5
            height[t, t] = t   
            Swalbe.fluid_dry!(fluid, dummy, height, t)
        end
        # println("This is fluid $fluid\nAnd sol $(reshape(sol2, 5, 25))")
        for t in 1:5
            @test all(fluid[t,:] .== sol2[t,:])
        end
    end

    @testset "Surface area liquid vapor and reduced energy" begin
        height = reshape(collect(1.0:25.0),5,5)
        dummy = ones(5,5,9)
        surface = ones(5,5)
        dx = ones(5,5)
        dy = ones(5,5)
        area_lv = [0.0]
        red_energy = [0.0]
        θ = fill(0.5, 5, 5)
        Swalbe.surfacearea!(area_lv, red_energy, height, θ, dx, dy, dummy, surface, 1)
        # Solutions gradient
        Xgrad = [-1.5 -1.5 -1.5 -1.5 -1.5;
                  1.0 1.0 1.0 1.0 1.0;
                  1.0 1.0 1.0 1.0 1.0;
                  1.0 1.0 1.0 1.0 1.0;
                 -1.5 -1.5 -1.5 -1.5 -1.5]

        Ygrad = [-7.5 5.0 5.0 5.0 -7.5;
                 -7.5 5.0 5.0 5.0 -7.5;
                 -7.5 5.0 5.0 5.0 -7.5;
                 -7.5 5.0 5.0 5.0 -7.5;
                 -7.5 5.0 5.0 5.0 -7.5]
        
        solsurf = 0.0
        solsurf = sum(sqrt.(Xgrad.^2 .+ Ygrad.^2 .+ 1))
        solener = 0.0
        solener = solsurf - sum(cospi.(θ))
        println("area is $(area_lv)")
        @test solsurf ≈ area_lv[1] atol=1e-10
        @test solener ≈ red_energy[1] atol=1e-10
        # With a non vanishing cosine
        θ = fill(1/9, 5, 5)
        Swalbe.surfacearea!(area_lv, red_energy, height, θ, dx, dy, dummy, surface, 1)
        solener = solsurf - sum(cospi.(θ))
        @test solsurf ≈ area_lv[1] atol=1e-10
        @test solener ≈ red_energy[1] atol=1e-10
        # More time steps
        d1 = zeros(5,5)
        d2 = zeros(5,5)
        area_lv = zeros(5)
        red_energy = zeros(5)
        surf_cal = zeros(5)
        ener_cal = zeros(5)
        for t in 1:5
            height[t,t] = t
            θ[t,t] = 0.0
            Swalbe.∇f_simple!(d1, d2, height, dummy)
            Swalbe.surfacearea!(area_lv, red_energy, height, θ, dx, dy, dummy, surface, t)
            surf_cal[t] = sum(sqrt.(d1.^2 .+ d2.^2 .+ 1))
            ener_cal[t] = surf_cal[t] - sum(cospi.(θ))
        end
        @test all(area_lv .== surf_cal)
        @test all(red_energy .== ener_cal)
    end
end