# =============================================================================
# Validation of SchlessingerEnsembleExtrapolator
# -----------------------------------------------------------------------------
# Ensemble of Schlessinger continued fractions, each built from the data with
# a distinct, seeded leave-out set. Evaluation: componentwise median, members
# farther than num_mad · MAD from it are dropped, the rest averaged; a warning
# is logged if the relative spread exceeds warn_tol.
#
# Run from the JuliaNumerics directory with
#   julia test/test_schlessinger_ensemble.jl [<path to Numerics/src>]
# =============================================================================

using Test, Logging

const NUMERICS_SRC = length(ARGS) >= 1 ? ARGS[1] : joinpath(@__DIR__, "..", "src")
include(joinpath(NUMERICS_SRC, "Numerics.jl"))

const SE = Numerics.Schlessinger


# Analytic, non-rational, real on the real axis; pole pair away from the data
const z0 = -1.5 + 1.0im
f_test(z) = log(z + 2.0) + 0.5 / ((z - z0) * (z - conj(z0)))

const n_data = 30
const x_data = collect(range(0.0, 3.0, length=n_data))
const f_data = real.(f_test.(x_data))

const num_leave_out = 2
const num_schlessingers = 40


@testset "SchlessingerEnsembleExtrapolator" begin

    ens = SE.SchlessingerEnsembleExtrapolator(num_schlessingers, num_leave_out, x_data, f_data; seed=42)

    @testset "members" begin
        @test length(ens.schl_array) == num_schlessingers

        # Each member keeps n - num_leave_out data points, ascending, with (x, f) pairs intact
        for schl in ens.schl_array
            @test length(schl.z_vals) == n_data - num_leave_out
            @test issorted(real.(schl.z_vals))
            @test all(x in x_data for x in schl.z_vals)
            @test schl.f_vals ≈ real.(f_test.(schl.z_vals)) rtol = 1e-14
        end

        # Leave-out sets are distinct
        @test allunique([Set(schl.z_vals) for schl in ens.schl_array])
    end

    @testset "seed" begin
        same = SE.SchlessingerEnsembleExtrapolator(num_schlessingers, num_leave_out, x_data, f_data; seed=42)
        other = SE.SchlessingerEnsembleExtrapolator(num_schlessingers, num_leave_out, x_data, f_data; seed=43)

        @test [s.z_vals for s in same.schl_array] == [s.z_vals for s in ens.schl_array]
        @test [s.z_vals for s in other.schl_array] != [s.z_vals for s in ens.schl_array]
    end

    @testset "too many members" begin
        @test_throws AssertionError SE.SchlessingerEnsembleExtrapolator(binomial(n_data, num_leave_out) + 1, num_leave_out, x_data, f_data)
    end

    @testset "reproduces analytic function" begin
        # Data nodes, between nodes, real extrapolation, complex plane
        for z in (x_data[7], 1.2345, 4.0, 1.5 + 0.8im, -0.5 - 0.5im)
            val = @test_logs min_level=Logging.Warn SE.value_at(z, ens)
            @test val ≈ f_test(z) rtol = 1e-6
        end
    end

    @testset "disagreement warning" begin
        # Noisy data: members disagree away from the data
        noisy = f_data .+ 1e-2 .* sin.(37.0 .* x_data)
        ens_noisy = SE.SchlessingerEnsembleExtrapolator(num_schlessingers, num_leave_out, x_data, noisy; seed=42)

        @test_logs (:warn, "Schlessinger ensemble disagrees") SE.value_at(6.0 + 2.0im, ens_noisy)
    end
end
