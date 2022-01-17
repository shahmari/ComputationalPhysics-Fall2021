
using Statistics: mean

mutable struct MDSystem{DT<:AbstractFloat}
    l::DT
    U::DT
    K::DT
    T::DT
    P::DT
    N::Integer
    r::Matrix{DT}
    v::Matrix{DT}
    function MDSystem(N::Integer, T₀::DT, l::DT) where {DT<:AbstractFloat}
        r = rand(DT, 2, N) * l
        v = rand(DT, 2, N) .- 0.5
        v .-= mean(v, dims = 2)
        v *= √(T₀ ÷ mean(v .^ 2))

        Σrterm = 0.0

        for i ∈ 1:(N-1)
            for j ∈ i+1:N
                Δr = abs(r[:, i] - r[:, j])
                if Δr > 2^(1 / 6)
                    Σrterm += ((Δr)^-12) - ((Δr)^-6)
                end
            end
        end

        U = 4 * Σrterm

        K = sum(v .^ 2) / 2

        T = K / N
        return new{DT}(l, U, K, T, P)
    end
end

function update_potential!(sys::MDSystem)
    Σrterm = 0.0

    for i ∈ 1:(sys.N-1)
        for j ∈ i+1:sys.N
            Δr = √sum(abs.(r[:, i] - r[:, j]) .^ 2)
            if Δr < 2^(1 / 6)
                Σrterm += ((Δr)^-12) - ((Δr)^-6)
            end
        end
    end

    U = 4 * Σrterm

    sys.U = U
end

function update_kinetic!(sys::MDSystem)
    sys.K = sum(sys.v .^ 2) / 2
end

function update_temperature!(sys::MDSystem)
    sys.T = sum(sys.v .^ 2) / (2 * sys.N)
end

function update_pressure!(sys::MDSystem)
    Σv² = sum(sys.v .^ 2)
    Σrterm = 0.0

    for i ∈ 1:(sys.N-1)
        for j ∈ i+1:sys.N
            Δr = √sum(abs.(r[:, i] - r[:, j]) .^ 2)
            if Δr < 2^(1 / 6)
                Σrterm += ((Δr)^-12) - 0.5 * ((Δr)^-6)
            end
        end
    end
    sys.P = ((48 * Σrterm) + Σv²) / (2 * sys.l * sys.l)
end

