
using Statistics: mean

mutable struct MDSystem{DT<:AbstractFloat}
    h::DT
    l::DT
    U::DT
    K::DT
    T::DT
    P::DT
    N::Integer
    f::Matrix{DT}
    r::Matrix{DT}
    v::Matrix{DT}
    function MDSystem(N::Integer, T₀::DT, l::DT, h) where {DT<:AbstractFloat}
        r = rand(DT, 2, N) * l
        v = rand(DT, 2, N) .- 0.5
        v .-= mean(v, dims = 2)
        v *= √(T₀ ÷ mean(v .^ 2))
        f = Matrix{DT}(undef, 2, N)
        return new{DT}(h, DT(NaN), DT(NaN), DT(NaN), DT(NaN), DT(NaN), N, f, r, v)
    end
end

function update_force!(sys::MDSystem)
    for i ∈ 1:sys.N
        fᵢ = sys.DT[0.0, 0.0]
        for j ∈ setdiff(1:sys.N, i)
            Δr = √sum(abs.(sys.r[:, i] - sys.r[:, j]) .^ 2)
            fᵢ += 48 * (((Δr)^-14) - 0.5 * ((Δr)^-8)) * abs.(sys.r[:, i] - sys.r[:, j])
        end
        sys.f[:, i] = fᵢ
    end
end

function get_force(sys::MDSystem, i)
    fᵢ = sys.DT[0.0, 0.0]
    for j ∈ setdiff(1:sys.N, i)
        Δr = √sum(abs.(sys.r[:, i] - sys.r[:, j]) .^ 2)
        fᵢ += 48 * (((Δr)^-14) - 0.5 * ((Δr)^-8)) * abs.(sys.r[:, i] - sys.r[:, j])
    end
    return fᵢ
end

function update_phase!(sys::MDSystem)
    for i ∈ 1:sys.N
        f₁ = sys.f[:, i]
        sys.r[:, i] += sys.h * sys.v[:, i] + 0.5 * sys.h * sys.h * f₁
        f₂ = get_force(sys, i)
        sys.f[:, i] = f₂
        sys.v[:, i] += sys.h * (f₁ + f₂) / 2
    end
end

function update_potential!(sys::MDSystem)
    Σrterm = 0.0

    for i ∈ 1:(sys.N-1)
        for j ∈ i+1:sys.N
            Δr = √sum(abs.(sys.r[:, i] - sys.r[:, j]) .^ 2)
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
            Δr = √sum(abs.(sys.r[:, i] - sys.r[:, j]) .^ 2)
            if Δr < 2^(1 / 6)
                Σrterm += ((Δr)^-12) - 0.5 * ((Δr)^-6)
            end
        end
    end
    sys.P = ((48 * Σrterm) + Σv²) / (2 * sys.l * sys.l)
end

MDSys = MDSystem(100, 200.0, 1.0)

MDSys.P

update_pressure!(MDSys)


