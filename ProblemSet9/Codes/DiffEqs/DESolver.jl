module DiffEqSolver

export EulerDES

function EulerDES(; ẋ::Function, x₀::Vector{Float64}, t₀::Float64, t₁::Float64, h::Float64)
    tcoll = collect(t₀:h:t₁)
    xcoll = Float64[x₀]
    vcoll = Float64[ẋ(t, x₀)]
    x = x₀
    for t ∈ tcoll
        x += h * ẋ(t, x)
        push!(xcoll, x)
        push!(vcoll, ẋ(t, x))
    end
    return tcoll, xcoll, vcoll
end

function EulerYaghoubDES(; ẋ::Function, x₀::Vector{Float64}, t₀::Float64, t₁::Float64, h::Float64)
    x = x₀
    v = v₀ = ẋ(t, x)
    tcoll = collect(t₀:h:t₁)
    xcoll = Float64[x₀]
    vcoll = Float64[v₀]
    for t ∈ tcoll
        x += h * v
        push!(xcoll, x)
        v = ẋ(t, x)
        push!(vcoll, v)
    end
    return tcoll, xcoll, vcoll
end

function EulerCromerDES(; ẍ::Function, x₀::Vector{Float64}, v₀::Vector{Float64}, t₀::Float64, t₁::Float64, h::Float64)
    tcoll = collect(t₀:h:t₁)
    xcoll = Float64[x₀]
    vcoll = Float64[v₀]
    x = x₀
    v = v₀
    for t ∈ tcoll
        v += h * ẍ(t, x)
        x += h * v
        push!(xcoll, x)
        push!(vcoll, v)
    end

    return ts, xs, vs
end

end