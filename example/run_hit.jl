using Revise
using HIT, WaterLily, CUDA, Printf, LaTeXStrings, Plots, Random
import WaterLily: dot, sgs!, size_u, @loop, @inside, inside, inside_u
Random.seed!(99) # seed random turbulence generator

function forced_hit_sgs!(flow, t; νₜ, S, C=0.16)
    # Forced injection
    forced_hit!(flow, t)
    # SGS modelling
    sgs!(flow, t; νₜ, S, C)
end
smagorinsky(I::CartesianIndex{m} where m, S; C=0.2) = @views C^2*sqrt(dot(S[I,:,:],S[I,:,:]))
function hit(N, M; cbc_path="data/cbc_spectrum.dat", ν=1e-6, L=1, U=1, mem=Array, T=Float32)
    sim = Simulation((N,N,N), (0,0,0), L; U, ν, T, mem, perdir=(1,2,3))
    u0, v0, w0 = generate_hit(L,N,M; cbc_path, mem)
    Ni,d = size_u(sim.flow.u)
    WaterLily.@loop sim.flow.u[I,1] = u0[I-CartesianIndex(2,1,1)] over I in WaterLily.inside_u(Ni,1)
    WaterLily.@loop sim.flow.u[I,2] = v0[I-CartesianIndex(1,2,1)] over I in WaterLily.inside_u(Ni,2)
    WaterLily.@loop sim.flow.u[I,3] = w0[I-CartesianIndex(1,1,2)] over I in WaterLily.inside_u(Ni,3)
    WaterLily.perBC!(sim.flow.u, sim.flow.perdir)
    return sim, stack([u0,v0,w0])
end
function run_hit(p, mem; udf=nothing, C=0.2, Re=1600, T=Float32, t_max=20.0, verbose=true)
    sim = hit(p, mem; Re, T)
    N,n = size_u(sim.flow.u)
    S = zeros(T, N..., n, n) |> mem # working array holding a tensor for each cell
    while WaterLily.time(sim.flow)*(sim.U/sim.L) < t_max
        sim_step!(sim; remeasure=false, udf, νₜ=smagorinsky, S, C)
        verbose && println("tU/L=", round(WaterLily.time(sim.flow)*(sim.U/sim.L), digits=4), ", Δt=",round(sim.flow.Δt[end], digits=3))

        kei, zi = kez!(sim.flow.σ, sim.flow.u, sim.flow.ν, sim.L)
        push!(ke, kei); push!(z, zi); push!(t, WaterLily.time(sim.flow)*(sim.U/sim.L))
    end
    return sim, ke, z, t
end

L = 9*2π/100 # reference length
N = 2^6 # cells per direction
M = 2^10 # number of modes for initial condition
ν = 1e-5
t0, t1, t2 = 42.0, 98.0, 171.0
cbc_length = 0.0508 # length scale related to grid size
cbc_velocity = 10.0 # velocity scale related to inflow
t_norm = cbc_length/cbc_velocity # convective time scale
T = Float32
mem = CuArray
C = 0.18 # Smagorinsky constant, C=0.18 typically
C_str = @sprintf("%2.2f", C)
udf = C > 0 ? sgs! : nothing

cbc_path = joinpath(string(@__DIR__), "../data/cbc_spectrum.dat")
set_plots_style!(; linewidth=2)

function main()
    sim, u = hit(N, M; cbc_path, ν, L, mem, T)
    t_str = @sprintf("%2.2f", t0)
    p = plot_spectra!(Plots.plot(), sim.flow.u|>Array, L, N;
        cbc_path, cbc_t=1, fig_path="plots/Ek_N$(N)_M$(M)_C$(C_str)_t$(t_str).pdf", label=L"t=%$t_str"
    )

    N1,n = size_u(sim.flow.u)
    S = zeros(T, N1..., n, n) |> mem # working array holding a tensor for each cell

    sim_step!(sim, t1-t0; verbose=true, remeasure=false, udf, νₜ=smagorinsky, S, C)
    t_str = @sprintf("%2.2f", sim_time(sim)+t0)
    p = plot_spectra!(p, sim.flow.u|>Array, L, N;
        cbc_path, cbc_t=2, fig_path="plots/Ek_N$(N)_M$(M)_C$(C_str)_t$t_str.pdf", label=L"t=%$t_str"
    )

    sim_step!(sim, sim_time(sim)+(t2-t1); verbose=true, remeasure=false, udf, νₜ=smagorinsky, S, C)
    t_str = @sprintf("%2.2f", sim_time(sim)+t0)
    p = plot_spectra!(p, sim.flow.u|>Array, L, N;
        cbc_path, cbc_t=3, fig_path="plots/Ek_N$(N)_M$(M)_C$(C_str)_t$t_str.pdf", label=L"t=%$t_str"
    )

    return sim, p
end

sim, p = main();