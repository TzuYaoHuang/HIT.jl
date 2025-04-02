using Revise
using HIT, Printf, LaTeXStrings, JLD2, FFTW, Plots
using WaterLily: Flow
using WaterLily
import EllipsisNotation: Ellipsis; const dots = Ellipsis()


M = 5.08/100 # grid size [m]
L = 9*2π/100 # length of HIT cube [m], L = 11M
velocity_scale = 10 # velocity related to the bulk flow (U₀ in paper) [m/s]

T = Float32 # run with single (Float32) or double (Float64) precision
mem = Array # run on CPU (Array) or GPU (CuArray)
N = 2^7 # cells per direction
cbc_path = joinpath(@__DIR__, "data", "cbc_spectrum.dat")

L_t = (L,L,L)
N_t = (N,N,N)
D = length(N_t)

sim = Simulation(N_t,(0,0,0),velocity_scale)

t_str = "42.00"
load!(sim.flow, joinpath(@__DIR__, "data/", "flow_N128_t$(t_str).jld2"))

# Filter
u = sim.flow.u[WaterLily.inside_u(sim.flow.u),:]
u_filtered = filter_sharp(u, N÷2) # filter at target grid size

f1, ax1 = σ_contour(u[:,:,N÷2,1])
f2, ax2 = σ_contour(u_filtered[:,:,N÷2,1])

set_plots_style!(; linewidth=2)
p = plot_spectra!(Plots.plot(dpi=600, title=L"N=%$(N)"), L, N, u;
    cbc_path, cbc_t=1, label=L"N=%$(N)"
)
p1 = plot_spectra!(p, L, N, u_filtered;
    cbc_path, cbc_t=1, label=L"N=%$(Nc)", fig_path="filtered.pdf"
)
