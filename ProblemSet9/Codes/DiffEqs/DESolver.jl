module DiffEqsSolver

export EulerDES, EulerYaghoubDES, EulerCromerDES, RK2DES, EulerMidPiontDES

function EulerDES(; ẋ::Function, x₀::Vector{Float64}, t₀::Float64, t₁::Float64, h::Float64)
    tcoll = collect(t₀:h:t₁)
    xcoll = [x₀]
    x = x₀
    for t ∈ tcoll[2:end]
        x += h * ẋ(t, x)
        push!(xcoll, x)
    end
    return tcoll, xcoll
end

function EulerYaghoubDES(; ẋ::Function, x₀::Vector{Float64}, t₀::Float64, t₁::Float64, h::Float64)
    tcoll = collect(t₀:h:t₁)
    xcoll = [x₀]
    x = x₀
    ẋval = ẋ(t, x₀)
    for t ∈ tcoll[2:end]
        x += h * ẋval
        ẋval = ẋ(t, x)
        push!(xcoll, x)
    end
    return tcoll, xcoll
end

function RK2DES(; ẋ::Function, x₀::Vector{Float64}, t₀::Float64, t₁::Float64, h::Float64)
    tcoll = collect(t₀:h:t₁)
    xcoll = [x₀]
    x = x₀
    for t ∈ tcoll[2:end]
        K₁ = ẋ(t, x) * h
        K₂ = h * ẋ(t + h / 2, x + (K₁ * h) / 2)
        x += K₂
        push!(xcoll, x)
    end
    return tcoll, xcoll
end

function EulerCromerDES(; ẍ::Function, x₀::Float64, t₀::Float64, t₁::Float64, h::Float64)
    tcoll = collect(t₀:h:t₁)
    xcoll = Float64[x₀]
    vcoll = Float64[v₀]
    x = x₀
    v = v₀
    for t ∈ tcoll[2:end]
        v += h * ẍ(t, x)
        x += h * v
        push!(xcoll, x)
        push!(vcoll, v)
    end

    return tcoll, xcoll, vcoll
end

function EulerMidPiontDES(; ẍ::Function, x₀::Float64, v₀::Float64, t₀::Float64, t₁::Float64, h::Float64)
    tcoll = collect(t₀:h:t₁)
    xcoll = Float64[x₀]
    vcoll = Float64[v₀-(ẍ(t₀, initx)*h/2)]
    x = x₀
    v = v₀

    for t ∈ tcoll[2:end]
        v += step * ẍ(t, x)
        x += step * v
    
        push!(xcoll, x)
        push!(vcoll, v)
    end

    return tcoll, xcoll, vcoll
end

end