using DFTK
using HDF5
using Libxc
using LinearAlgebra

if !isfile("../sigma.hdf5")
    lattice = diagm([3.0, 0.0, 0.0])
    Si = ElementPsp(:Si, psp=load_psp("hgh/lda/si-q4"))
    model = model_PBE(lattice, [Si => [ones(3)/8, -ones(3)/8]],
                      spin_polarization=:collinear, temperature=0.01)

    Ecut = 2
    basis = PlaneWaveBasis(model, Ecut; kgrid=[1, 1, 1])

    ρ = guess_density(basis).real
    ρspin = rand(basis.fft_size...)
    ρspin[abs.(ρspin) .> ρ] .= 0.0
    ρspin = from_real(basis, ρspin)
    ρ     = from_real(basis, ρ)

    density = DFTK.LibxcDensity(basis, 1, ρ, ρspin)
    h5write("../sigma.hdf5", "ρ", density.ρ_real)
    h5write("../sigma.hdf5", "σ", density.σ_real)
end
ρ_real = reshape(h5read("../sigma.hdf5", "ρ"), :)
σ_real = reshape(h5read("../sigma.hdf5", "σ"), :)

@show ρ_real
@show σ_real
xc = Functional(:gga_c_pbe; n_spin=2)
zk = evaluate(xc; rho=ρ_real, sigma=σ_real).zk
@show zk
nothing
