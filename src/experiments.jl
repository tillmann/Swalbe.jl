"""
    dropletrelax()

Simulates an out of equilibrium droplet
"""
function dropletrelax(
    sys::SysConst, 
    device::String; 
    radius=20, 
    θ₀=1/6, 
    center=(sys.Lx÷2, sys.Ly÷2), 
    verbos=true, 
    T=Float64
)
    println("Simulating an out of equilibrium droplet")
    area = []
    fout, ftemp, feq, height, velx, vely, vsq, pressure, Fx, Fy, slipx, slipy, h∇px, h∇py = Swalbe.Sys(sys, device, false, T)
    Swalbe.singledroplet(height, radius, θ₀, center)
    Swalbe.equilibrium!(feq, height, velx, vely, vsq)
    ftemp .= feq
    for t in 1:sys.Tmax
        if t % sys.tdump == 0
            mass = 0.0
            mass = sum(height)
            difference = maximum(height) - minimum(height)
            if verbos
                println("Time step $t mass is $(round(mass, digits=3))")
            end
        end
        push!(area, length(findall(height .> 0.055)))
        Swalbe.filmpressure!(pressure, height, sys.γ, 1/9, sys.n, sys.m, sys.hmin, sys.hcrit)
        Swalbe.∇f!(h∇px, h∇py, pressure, height)
        Swalbe.slippage!(slipx, slipy, height, velx, vely, sys.δ, sys.μ)
        Fx .= h∇px .+ slipx
        Fy .= h∇py .+ slipy
        Swalbe.equilibrium!(feq, height, velx, vely, vsq)
        Swalbe.BGKandStream!(fout, feq, ftemp, -Fx, -Fy)
        Swalbe.moments!(height, velx, vely, fout)
    end
    return height, area
end