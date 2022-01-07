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

function VerletDES(; ẍ::Function, x₀::Float64, v₀::Float64, t₀::Float64, t₁::Float64, h::Float64)
    tcoll = collect(t₀:h:t₁)
    xcoll = Float64[x₀, x₀+v₀*h+ẍ(t₀, x₀)*(h^2)/2]
    vcoll = Float64[v₀, v₀+ẍ(t₀, x₀)*h]

    for i ∈ 3:length(tcoll)
        xᵢ = 2 * xcoll[i-1] - xcoll[i-2] + ẍ(tcoll[i-1], xcoll[i-1]) * (h^2)
        vᵢ = (xcoll[i] - xcoll[i-1]) / h
        push!(xcoll, xᵢ)
        push!(vcoll, vᵢ)
    end

    return tcoll, xcoll, vcoll
end

function VelVerDES(; ẍ::Function, x₀::Float64, v₀::Float64, t₀::Float64, t₁::Float64, h::Float64)
    tcoll = collect(t₀:h:t₁)
    xcoll = Float64[x₀]
    vcoll = Float64[v₀]
    x = x₀
    v = v₀

    for t ∈ tcoll[2:end]
        ẍₙ = ẍ(t, x)
        x += v * h + ẍₙ * (h^2) / 2
        ẍₙ₊₁ = ẍ(t + h, x)
        v += (ẍₙ + ẍₙ₊₁) * h / 2
        push!(xcoll, x)
        push!(vcoll, v)
    end

    return tcoll, xcoll, vcoll
end

end