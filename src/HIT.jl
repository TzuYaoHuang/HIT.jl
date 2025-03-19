module HIT
using DelimitedFiles, Interpolations, Distributions, Random, OutMacro, FFTW, CUDA
import EllipsisNotation: Ellipsis; const dots = Ellipsis()
import WaterLily: dot

include("util.jl")

export generate_hit, spectrum, cbc_spectrum, plot_spectra!, ω_viz, set_plots_style!

"""
    generate_hit(L,N,M; mem=Array)

Homogeneous isotropic turbulence initial condition.
    `L`: domain size
    `N`: number of grid points per direction
    `M`: number of modes
"""
function generate_hit(L,N,M; cbc_path="data/cbc_spectrum.dat", mem=Array)
    dx = L/N # cell size
    wn1 = 2π/L # smallest wavenumber represented by this spectrum, determined here from cbc spectrum properties

    # Compute random angles in unit sphere
    ν = rand(Uniform(0,1), M)
    ϕm = rand(Uniform(0,2π), M)
    θm = @. acos(2ν-1)
    ψm = rand(Uniform(-π/2,π/2), M)
    # Highest wave number that can be represented on this grid (nyquist limit)
    wnn = 2π/dx
    # Wavenumber step
    dk = (wnn - wn1) / M
    # Wavenumber at cell centers
    wn = @. wn1 + 0.5 * dk + $(collect(0:M-1)) * dk
    # Wavenumber vector from random angles
    kx = @. sin(θm) * cos(ϕm) * wn
    ky = @. sin(θm) * sin(ϕm) * wn
    kz = @. cos(θm) * wn
    # Create divergence vector
    ktx = @. sin(kx*dx/2)/dx
    kty = @. sin(ky*dx/2)/dx
    ktz = @. sin(kz*dx/2)/dx
    # Enforce Mass Conservation
    ϕm1 = rand(Uniform(0,2π), M)
    ν1 = rand(Uniform(0,1), M)
    θm1 = @. acos(2ν1-1)
    ζx = @. sin(θm1) * cos(ϕm1)
    ζy = @. sin(θm1) * sin(ϕm1)
    ζz = @. cos(θm1)
    sxm = @. ζy * ktz - ζz * kty
    sym = @. -(ζx * ktz - ζz * ktx)
    szm = @. ζx * kty - ζy * ktx
    smag = @. sqrt(sxm * sxm + sym * sym + szm * szm)
    @. sxm = sxm / smag
    @. sym = sym / smag
    @. szm = szm / smag

    # Verify that the wave vector and sigma are perpendicular
    # @assert isapprox(sum(dot(ktx, sxm) + dot(kty, sym) + dot(ktz, szm)), 0; atol=100eps(T)) "wave vector and sigma are not perpendicular"

    # Get CBC spectrum
    k_cbc, E = cbc_spectrum(cbc_path)
    # Generate turbulence at cell centers
    um = @. sqrt(E(wn)*dk)
    u,v,w = zeros(N,N,N), zeros(N,N,N), zeros(N,N,N)
    c = dx/2 .+ collect(0:N-1)*dx # cell centers

    arg, bmx, bmy, bmz = zeros(M), zeros(M), zeros(M), zeros(M)
    for k=1:N, j=1:N, i=1:N
        @. arg = kx * c[i] + ky * c[j] + kz * c[k] - ψm
        @. bmx = 2.0 * um * cos(arg - kx * dx / 2.0)
        @. bmy = 2.0 * um * cos(arg - ky * dx / 2.0)
        @. bmz = 2.0 * um * cos(arg - kz * dx / 2.0)

        u[i,j,k] = dot(bmx, sxm)
        v[i,j,k] = dot(bmy, sym)
        w[i,j,k] = dot(bmz, szm)
    end
    return (u, v, w) .|> mem
end

spectrum(u, L::Number) = spectrum(u, Tuple(L for i in 1:last(size(u))))
function spectrum(u, L::Tuple)
    N,d = size(u)[1:end-1], last(size(u))
    @assert length(L) == d

    dx = L ./ N
    k0 = 2π ./ L
    knorm = mean(k0)
    wn_vec_out = collect(knorm * i for i in 0:N[1]-1)
    wn_vec = collect(fftfreq(N[i], dx[i]) * N[i]/dx[i] for i in 1:d) # or wave numbers, 0..(N-2)/2,-N/2,...,-1, vcat(0:(N[i]-2)/2,-N[i]/2:-1)
    r_wn = collect(sqrt(sum(wn_vec[i][I[i]]^2 for i in 1:d)) for I in CartesianIndices(N))
    for i in 1:d
        u[dots,i] .-= mean(u[dots,i])
    end
    uk = collect(fft(u[dots,i])/prod(N) for i in 1:d)
    tke = 0.5sum(uk[dots,i].*conj(uk[dots,i]) for i in 1:d) |> real # TKE distributed across n-dimensional wavenumber space
    tke_sum = zeros(length(wn_vec_out)) # spherically integrated TKE with M modes of resolution

    for ijk in eachindex(tke)
        k = round(Int, r_wn[ijk])
        k < 1 && continue
        tke_sum[k] += tke[ijk]
    end
    return wn_vec_out, tke_sum./knorm
end

end # module HIT
