using Revise
using HIT, WaterLily, CUDA, Printf, LaTeXStrings, Plots, Random
import WaterLily: dot, sgs!, size_u, @loop, @inside, inside, inside_u, CFL
Random.seed!(99) # seed random turbulence generator


smagorinsky(I::CartesianIndex{m} where m, S; C=0.2) = @views C^2*sqrt(dot(S[I,:,:],S[I,:,:]))

function hit(L, N, M; length_scale=1, velocity_scale=1, cbc_path="data/cbc_spectrum.dat", ν=1e-6, mem=Array, T=Float32)
    sim = Simulation((N,N,N), (0,0,0), length_scale; U=velocity_scale, ν, T, mem, perdir=(1,2,3))
    u0, v0, w0 = generate_hit(L,N,M; cbc_path, mem)
    Ni,d = size_u(sim.flow.u)
    WaterLily.@loop sim.flow.u[I,1] = u0[I-CartesianIndex(2,1,1)] over I in WaterLily.inside_u(Ni,1)
    WaterLily.@loop sim.flow.u[I,2] = v0[I-CartesianIndex(1,2,1)] over I in WaterLily.inside_u(Ni,2)
    WaterLily.@loop sim.flow.u[I,3] = w0[I-CartesianIndex(1,1,2)] over I in WaterLily.inside_u(Ni,3)
    WaterLily.perBC!(sim.flow.u, sim.flow.perdir)
    return sim
end

M = 5.08/100 # grid size [m]
L = 9*2π/100 # length of HIT cube # L ≈ 11grid_size
N = 2^7 # cells per direction
modes = 2^10 # number of modes for initial condition
ν = 1.5e-5 # same as Bae et al
length_scale = N / (L/M) # for CTU use M, which is resolved with N/11 (11 Ms per N)
velocity_scale = 10 # 10.0 velocity scale related to inflow (Uo in paper)
t0_ctu, t1_ctu, t2_ctu = 42.0, 98.0, 171.0 # in ctu, t_ctu=length_scale/velocity_scale = M/U
T = Float32
mem = CuArray
C = 0.18 # Smagorinsky constant, C=0.18 typically
C_str = @sprintf("%2.2f", C)
udf = C > 0 ? sgs! : nothing

cbc_path = joinpath(string(@__DIR__), "../data/cbc_spectrum.dat")
set_plots_style!(; linewidth=2)

function main()
    println("Running with SGS: $udf")
    sim = hit(L, N, modes; length_scale, velocity_scale, ν, mem, T)
    u_inside = @views sim.flow.u[inside_u(sim.flow.u),:]
    t_str = @sprintf("%2.2f", t0_ctu)
    p = plot_spectra!(Plots.plot(), L, N, u_inside|>Array;
        cbc_path, cbc_t=1, label=L"t=%$t_str"
    )

    N1,n = size_u(sim.flow.u)
    S = zeros(T, N1..., n, n) |> mem # working array holding a tensor for each cell

    sim_step!(sim, t1_ctu-t0_ctu; verbose=true, remeasure=false, udf, νₜ=smagorinsky, S, C)
    t_str = @sprintf("%2.2f", sim_time(sim)+t0_ctu)
    p = plot_spectra!(p, L, N, u_inside|>Array;
        cbc_path, cbc_t=2, label=L"t=%$t_str"
    )

    sim_step!(sim, sim_time(sim)+(t2_ctu-t1_ctu); verbose=true, remeasure=false, udf, νₜ=smagorinsky, S, C)
    t_str = @sprintf("%2.2f", sim_time(sim)+t0_ctu)
    p = plot_spectra!(p, L, N, u_inside|>Array;
        cbc_path, cbc_t=3, fig_path="plots/Ek_N$(N)_modes$(modes)_C$(C_str)_t$t_str.pdf", label=L"t=%$t_str"
    )

    return sim, p
end

sim, p = main();