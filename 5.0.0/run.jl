using Test
using DFTK

@testset "Spin-broken silicon setup relaxes to spin-paired ground state" begin
    function run_silicon(spin_polarization)
        lattice = [0.0  5.131570667152971 5.131570667152971;
                   5.131570667152971 0.0 5.131570667152971;
                   5.131570667152971 5.131570667152971  0.0]
        Si = ElementPsp(:Si, psp=load_psp("hgh/lda/si-q4"))
        model = model_PBE(lattice, [Si => [ones(3)/8, -ones(3)/8]],
                          spin_polarization=spin_polarization, temperature=0.01)

        Ecut = 7
        basis = PlaneWaveBasis(model, Ecut; kgrid=[2, 2, 2])

        if spin_polarization == :collinear
            ρspin = from_real(basis, 1.0rand(basis.fft_size...))
        else
            ρspin = nothing
        end
        self_consistent_field(basis, tol=5e-6, ρspin=ρspin, n_bands=10);
    end

    scfres        = run_silicon(:none)
    scfres_broken = run_silicon(:collinear)
    εbroken       = scfres_broken.eigenvalues

    @test scfres.energies.total ≈ scfres_broken.energies.total atol=1e-5
    absmax(x) = maximum(abs, x)
    for (ik, kpt) in enumerate(scfres.basis.kpoints)
        kequiv = findall(kbr -> kbr.coordinate == kpt.coordinate, scfres_broken.basis.kpoints)

        for ikb in kequiv
            @test scfres.eigenvalues[ik][1:10] ≈ εbroken[ikb][1:10] atol=5e-4 norm=absmax
        end
    end
end
