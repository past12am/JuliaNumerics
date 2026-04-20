module Schlessinger
    
    struct SchlessingerExtrapolator{Ta<:Number}
        n::Int

        z_vals::AbstractArray{<:Number}
        f_vals::AbstractArray{<:Number}

        a::AbstractArray{Ta}

        SchlessingerExtrapolator{Ta}(z_vals::AbstractArray{<:Number}, f_vals::AbstractArray{<:Number}) where {Ta<:Number} = begin
            n = length(z_vals)

            # Build reciprocal difference table
            zc = Float64.(z_vals)
            fc = Float64.(f_vals)

            rd = zeros(ComplexF64, n, n)
            rd[1, :] .= fc

            for p in 1:n-1
                for j in p+1:n
                    rd[p+1, j] = (rd[p, p] / rd[p, j] - 1) / (zc[j] - zc[p])
                end
            end

            a = [rd[p+1, p+1] for p in 1:n-1]

            return new{Ta}(n, z_vals, f_vals, a)
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


end