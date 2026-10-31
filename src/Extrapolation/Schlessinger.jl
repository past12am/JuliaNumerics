module Schlessinger
    
    struct SchlessingerExtrapolator{Ta<:Number}
        n::Int

        z_vals::Vector{Ta}
        f_vals::Vector{Ta}

        a::Vector{Ta}

        SchlessingerExtrapolator{Ta}(z_vals::AbstractArray{Ta}, f_vals::AbstractArray{Ta}) where {Ta<:Number} = begin
            n = length(z_vals)

            # Build reciprocal difference table
            zc = Ta.(z_vals)
            fc = Ta.(f_vals)

            rd = zeros(ComplexF64, n, n)
            rd[1, :] .= fc

            for p in 1:n-1
                for j in p+1:n
                    rd[p+1, j] = (rd[p, p] / rd[p, j] - 1) / (zc[j] - zc[p])
                end
            end

            a = [rd[p+1, p+1] for p in 1:n-1]

            return new{Ta}(n, zc, fc, a)
        end
    end


    function value_at(z::Number, schlessinger::SchlessingerExtrapolator{T}) where {T<:Number}
        zc = ComplexF64(z)

        den::ComplexF64 = 1.0 + 0.0im
        for k in length(schlessinger.a):-1:1
            den = 1.0 + schlessinger.a[k] * (zc - schlessinger.z_vals[k]) / den
        end

        return schlessinger.f_vals[1] / den
    end


    using Random, Combinatorics, Statistics

    struct SchlessingerEnsembleExtrapolator{Ta<:Number}

        num_schlessingers::Int
        num_leave_out::Int

        schl_array::Vector{SchlessingerExtrapolator{Ta}}

        SchlessingerEnsembleExtrapolator(num_schlessingers::Int, num_leave_out::Int, x_vals::Vector{Ta}, f_vals::Vector{Ta}; seed=1234) where {Ta<:Number} = begin
            
            n = length(f_vals)
            @assert num_schlessingers <= binomial(n, num_leave_out)
            
            # Distinct leave-out sets, seeded order
            leave_outs = shuffle(Xoshiro(seed), collect(combinations(1:n, num_leave_out)))[1:num_schlessingers]

            schl_array = Vector{SchlessingerExtrapolator{Ta}}(undef, num_schlessingers)
            for (i, left_out) in enumerate(leave_outs)
                idx = setdiff(1:n, left_out)
                schl_array[i] = SchlessingerExtrapolator{Ta}(x_vals[idx], f_vals[idx])
            end

            new{Ta}(num_schlessingers, num_leave_out, schl_array)
        end
    end

    function value_at(z::Number, ensemble::SchlessingerEnsembleExtrapolator; num_mad=3.0, warn_tol=1e-2, do_warn=false)
        vals = [value_at(z, schl) for schl in ensemble.schl_array]

        # Componentwise median and spread
        med = median(real.(vals)) + im * median(imag.(vals))

        dist = abs.(vals .- med)
        mad = median(dist)

        # Warn on member disagreement
        if do_warn && mad > warn_tol * abs(med)
            @warn "Schlessinger ensemble disagrees" z rel_spread=mad/abs(med) num_kept=count(dist .<= num_mad * mad)
        end

        # Mean of members within num_mad · MAD
        return mean(vals[dist .<= num_mad * mad])
    end
end