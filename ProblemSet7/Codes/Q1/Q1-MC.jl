module MCIntegrate

export SMC, IMC

using Distributions, Statistics

function SMC(; f::Function, X₀::Real, X₁::Real, Ŋ::Integer = 100000)
    ΔX = X₁ - X₀
    fₓ = f.(ΔX .* rand(Ŋ) .+ X₀)
    Avg = mean(fₓ)
    ∫fₓdx = Avg * ΔX
    Δ = std(fₓ) / √Ŋ
    return ∫fₓdx, Δ
end

function IMC(; f::Function, g::Function, ∫gₓdx::Real, Đₓ::Function, Ŋ::Integer = 100000)
    fₓ╱gₓ = x -> f(x) / g(x)
    Samples = fₓ╱gₓ.(Đₓ.(rand(Ŋ)))
    Avg = mean(Samples)
    σ = std(Samples)
    ∫fₓdx = Avg * ∫gₓdx
    Δ = ∫gₓdx * σ / √Ŋ
    return ∫fₓdx, Δ
end

end
