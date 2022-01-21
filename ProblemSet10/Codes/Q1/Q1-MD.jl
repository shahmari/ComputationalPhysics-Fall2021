
using Statistics: mean
using Plots

mutable struct MDSystem{DT<:AbstractFloat}
    l::DT
    h::DT
    U::DT
    K::DT
    T::DT
    P::DT
    N::Integer
    f::Matrix{DT}
    r::Matrix{DT}
    v::Matrix{DT}
    function MDSystem(N::Integer, T₀::DT, h::DT, l::DT) where {DT<:AbstractFloat}
        r = rand(DT, 2, N) * l
        v = rand(DT, 2, N) .- 0.5
        v .-= mean(v, dims = 2)
        v *= √(T₀ ÷ mean(v .^ 2))
        f = Matrix{DT}(undef, 2, N)
        return new{DT}(l, h, DT(NaN), DT(NaN), DT(NaN), DT(NaN), N, f, r, v)
    end
end

function init(; N::Integer, T₀::DT, h::DT, l::DT) where {DT<:AbstractFloat}
    sys = MDSystem(N, T₀, h, l)
    alignleft!(sys)
    update_forces!(sys)
    update_phase!(sys)
    check_boundary!(sys)
    update_potential!(sys)
    update_kinetic!(sys)
    update_temperature!(sys)
    update_pressure!(sys)
    return sys
end

function simulate!(sys::MDSystem)
    update_phase!(sys)
    check_boundary!(sys)
    update_potential!(sys)
    update_kinetic!(sys)
    update_temperature!(sys)
    update_pressure!(sys)
end

function alignleft!(sys::MDSystem)    #Thanks to my dearest friend for his beautiful method in setting the initial conditions SLH Hero
    xcount = ceil(Integer, (2 * sys.N)^(1 / 2) / 2)
    ycount = 2 * xcount

    xs = collect(range(0.01 * sys.l / 2, 0.99 * sys.l / 2, length = xcount))
    sys.r[1, :] .= repeat(xs, outer = ycount)[1:sys.N]

    ys = collect(range(0.01 * sys.l, 0.99 * sys.l, length = ycount))
    sys.r[2, :] .= repeat(ys, inner = xcount, outer = 1)[1:sys.N]
end

function update_forces!(sys::MDSystem)
    for i ∈ 1:sys.N
        sys.f[:, i] *= 0
        for j ∈ setdiff(1:sys.N, i)
            rmirror = sys.r[:, j]
            if sys.r[1, i] - rmirror[1] > sys.l / 2
                rmirror[1] -= sys.l
            elseif sys.r[1, i] - rmirror[1] < -sys.l / 2
                rmirror[1] += sys.l
            end
            if sys.r[2, i] - rmirror[2] > sys.l / 2
                rmirror[2] -= sys.l
            elseif sys.r[2, i] - rmirror[2] < -sys.l / 2
                rmirror[2] += sys.l
            end
            Δr = √sum((sys.r[:, i] - rmirror) .^ 2)
            sys.f[:, i] += 48 * (((Δr)^-14) - 0.5 * ((Δr)^-8)) * (sys.r[:, i] - rmirror)
        end
    end
end


function update_force!(sys::MDSystem, i::Integer)
    sys.f[:, i] *= 0
    for j ∈ setdiff(1:sys.N, i)
        rmirror = sys.r[:, j]
        if sys.r[1, i] - rmirror[1] > sys.l / 2
            rmirror[1] -= sys.l
        elseif sys.r[1, i] - rmirror[1] < -sys.l / 2
            rmirror[1] += sys.l
        end
        if sys.r[2, i] - rmirror[2] > sys.l / 2
            rmirror[2] -= sys.l
        elseif sys.r[2, i] - rmirror[2] < -sys.l / 2
            rmirror[2] += sys.l
        end
        Δr = √sum((sys.r[:, i] - rmirror) .^ 2)
        sys.f[:, i] += 48 * (((Δr)^-14) - 0.5 * ((Δr)^-8)) * (sys.r[:, i] - rmirror)
    end
end

"""
important note:
    update_forces!() will update the force of the whole system,
    but update_force! only updates the selected particle.
"""

function update_phase!(sys::MDSystem)
    for i ∈ 1:sys.N
        f₁ = sys.f[:, i]
        sys.r[:, i] += sys.h * sys.v[:, i] + 0.5 * sys.h * sys.h * f₁
        update_force!(sys, i)
        f₂ = sys.f[:, i]
        sys.v[:, i] += sys.h * (f₁ + f₂) / 2
    end
end

function check_boundary!(sys::MDSystem)
    sys.r = (sys.r .+ sys.l) .% sys.l
end

function update_potential!(sys::MDSystem)
    Σrterm = 0.0

    for i ∈ 1:(sys.N-1)
        for j ∈ i+1:sys.N
            rmirror = sys.r[:, j]
            if sys.r[1, i] - rmirror[1] > 0.5
                rmirror[1] -= 1
            elseif sys.r[1, i] - rmirror[1] < -0.5
                rmirror[1] += 1
            end
            if sys.r[2, i] - rmirror[2] > 0.5
                rmirror[2] -= 1
            elseif sys.r[2, i] - rmirror[2] < -0.5
                rmirror[2] += 1
            end
            Δr = √sum((sys.r[:, i] - rmirror) .^ 2)
            Σrterm += ((Δr)^-12) - ((Δr)^-6)
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
            if Δr < 0.5
                Σrterm += ((Δr)^-12) - 0.5 * ((Δr)^-6)
            end
        end
    end
    sys.P = ((48 * Σrterm) + Σv²) / (2 * sys.l * sys.l)
end



Parameters = Dict(:N => 50, :T₀ => 1.0, :h => 0.005, :l => 30.0)
sys = init(; Parameters...)

@gif for i ∈ 1:100
    simulate!(sys)
    scatter([Tuple(sys.r[:, n]) for n ∈ 1:sys.N], xlims = (-1, 30.0), ylims = (-1, 30.0))
end


scatter([Tuple(sys.r[:, n]) for n ∈ 1:sys.N], xlims = (-1, 30.0), ylims = (-1, 30.0))