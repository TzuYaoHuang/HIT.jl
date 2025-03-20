using Plots, StatsPlots, LaTeXStrings, CategoricalArrays, Printf, ColorSchemes
using GLMakie; GLMakie.activate!(inline=false)
using WaterLily

set_plots_style!(; fontsize=14, linewidth=1) = Plots.default(
    fontfamily = "Computer Modern",
    linewidth = linewidth,
    framestyle = :box,
    grid = true,
    minorgrid = true,
    left_margin = Plots.Measures.Length(:mm, 2),
    right_margin = Plots.Measures.Length(:mm, 2),
    bottom_margin = Plots.Measures.Length(:mm, 2),
    top_margin = Plots.Measures.Length(:mm, 2),
    titlefontsize = fontsize,
    legendfontsize = fontsize,
    tickfontsize = fontsize,
    labelfontsize = fontsize,
)

function cbc_spectrum(cbc_path="data/cbc_spectrum.dat", cbc_t=1)
    cbc_spectrum = readdlm(cbc_path)
    k_cbc = 100 * cbc_spectrum[:, 1]
    e_cbc = 1e-6 * cbc_spectrum[:, 1+cbc_t]
    E = interpolate((k_cbc,), e_cbc, Gridded(Linear()))
    return k_cbc, E
end

function plot_spectra!(p, L, N, u; cbc_path="cbc_spectrum.dat", cbc_t=1, fig_path="Ek.pdf", label=L"~%$(N[1])^3")
    if !isnothing(cbc_path)
        label_cbc = any(p.series_list[i][:label]=="CBC" for i in 1:length(p.series_list)) ? nothing : "CBC"
        k_cbc, E = cbc_spectrum(cbc_path, cbc_t)
        Plots.plot!(p, k_cbc, E(k_cbc), label=label_cbc, color=:black)
    end

    k, tke = spectrum(u, L)
    Plots.plot!(p, k[2:end], tke[2:end], label=label, marker=:circle, markersize=2, markevery=1,)
    Plots.vline!(p, [2π/(L/(N/2))], label=:none, ls=:dash, color=:purple)

    Plots.plot!(p, xaxis=:log10, yaxis=:log10, xlims=(10,1e3), ylims=(1e-6,1e-3),
        xlabel=L"\kappa", ylabel=L"E(\kappa)",framestyle=:box, grid=true, minorgrid=true,
        left_margin=Plots.Measures.Length(:mm, 0), bottom_margin=Plots.Measures.Length(:mm, 0), size=(900,600)
    )
    savefig(p, fig_path)
    println("Figure stored in $(fig_path)")
    return p
end

function ω!(cpu_array, sim)
    a,dt = sim.flow.σ,sim.L/sim.U
    WaterLily.@inside a[I] = WaterLily.ω_mag(I,sim.flow.u)
    copyto!(cpu_array, a[inside(a)]) # copy to CPU
end

function ω_viz(sim; t_end=nothing, dt=0.0025, video=false)
    function viz_step!(sim; dt)
        sim_step!(sim, sim_time(sim)+dt; remeasure=false, verbose=true)
        ω[] = ω!(dat,sim)
    end

    dat = sim.flow.σ[inside(sim.flow.σ)] |> Array; # CPU buffer array
    ω = ω!(dat, sim) |> Observable
    f = Figure(size=(1200,1200))
    ax = Axis3(f[1, 1]; aspect=:equal, limits=(1,N,1,N,1,N))
    hidedecorations!(ax)

    colormap = to_colormap(:plasma)
    colormap[1] = RGBAf(0,0,0,0)
    # volume!(ax, ω,  algorithm = :absorption, absorption=1f0, colormap=colormap)
    volume!(ax, ω, algorithm=:iso, colormap=:rainbow, isovalue=0.1)
    display(f)

    if !isnothing(t_end) # time loop for animation
        if video
            GLMakie.record(f, "hit.mp4", 1:Int((t_end-sim_time(sim))/dt); framerate=60, compression=5) do frame
                viz_step!(sim; dt)
            end
        else
            while sim_time(sim) < t_end
                viz_step!(sim; dt)
            end
        end
    end
end