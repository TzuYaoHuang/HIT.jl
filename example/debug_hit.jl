using Revise
using HIT, Printf, LaTeXStrings, Plots, Random
Random.seed!(99) # seed random turbulence generator

L = 9*2π/100 # reference length
N = 2^6 # cells per direction
M = 2^6 # number of modes for initial condition
ν = 1e-5
t0, t1, t2 = 42.0, 98.0, 171.0
cbc_length = 0.0508 # length scale related to grid size
cbc_velocity = 10.0 # velocity scale related to inflow
T = Float64
mem = Array

cbc_path = joinpath(string(@__DIR__), "../data/cbc_spectrum.dat")
set_plots_style!(; linewidth=2)

function main()
    u,v,w = generate_hit(L,N,M; cbc_path, mem)
    ui = stack([u,v,w])
    t_str = @sprintf("%2.2f", t0)
    p = plot_spectra!(Plots.plot(), L, N, ui;
        cbc_path, cbc_t=1, fig_path="plots/test3.pdf", label=L"t=%$t_str"
    )
    return u,v,w
end

u,v,w = main();