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


    # Allocation free value_at

    struct EnsembleScratch
        vals::Vector{ComplexF64}
        work::Vector{Float64}

        EnsembleScratch(max_members::Int) = new(zeros(ComplexF64, max_members), zeros(Float64, max_members))
    end

    # Median of x, reordering x; allocation-free (insertion sort, the ensembles are small). Same value as Statistics.median!
    function median_inplace!(x::AbstractVector{Float64})
        n::Int = length(x)
        for i = 1:n
            isnan(x[i]) && return NaN
        end

        for i = 2:n
            xi::Float64 = x[i]
            j::Int = i - 1
            while j >= 1 && x[j] > xi
                x[j + 1] = x[j]
                j -= 1
            end
            x[j + 1] = xi
        end

        m::Int = (n + 1) ÷ 2
        return isodd(n) ? x[m] : middle(x[m], x[m + 1])
    end

    function value_at!(scratch::EnsembleScratch, z::Number, ensemble::SchlessingerEnsembleExtrapolator; num_mad=3.0, warn_tol=1e-2, do_warn=false)
        
        # Init arrays as views on scratch
        n::Int = ensemble.num_schlessingers
        @assert n <= length(scratch.vals) "EnsembleScratch holds $(length(scratch.vals)) members, the ensemble has $n"

        vals = view(scratch.vals, 1:n)
        work = view(scratch.work, 1:n)

        
        # Set member values
        for i = 1:n
            vals[i] = value_at(z, ensemble.schl_array[i])
        end


        # Componentwise median and spread (work is reordered, vals not)
        for i = 1:n
            work[i] = real(vals[i])
        end
        med_re::Float64 = median_inplace!(work)

        for i = 1:n
            work[i] = imag(vals[i])
        end
        med_im::Float64 = median_inplace!(work)

        med::ComplexF64 = med_re + im * med_im

        for i = 1:n
            work[i] = abs(vals[i] - med)
        end
        mad::Float64 = median_inplace!(work)


        # Members within num_mad · MAD, moved to the front of vals (order kept)
        num_kept::Int = 0
        for i = 1:n
            if abs(vals[i] - med) <= num_mad * mad
                num_kept += 1
                vals[num_kept] = vals[i]
            end
        end
        

        # Warn on member disagreement
        if do_warn && mad > warn_tol * abs(med)
            @warn "Schlessinger ensemble disagrees" z rel_spread=mad/abs(med) num_kept=num_kept
        end

        
        # Mean of the kept members
        return mean(view(vals, 1:num_kept))
    end
end