using Plots, StatsBase, LaTeXStrings, JLD


function SMC(;f::Function, X₀::Integer, X₁::Integer, Ŋ::Integer = 100000)
    ΔX = X₁ - X₀
    fₓ = f.(ΔX.*rand(Ŋ) .+ X₀)
    Avg = mean(fₓ)
    ∫fₓdx = Avg * ΔX
    Δ = std(fₓ) / √Ŋ
    return ∫fₓdx, Δ
end

function IMC(;f::Function, g::Function, ∫gₓdx::Real, Đₓ::Function, Ŋ::Integer = 100000)
    fₓ╱gₓ = x -> f(x)/g(x)
    Samples = fₓ╱gₓ.(Đₓ.(rand(Ŋ)))
    Avg = mean(Samples)
    σ = std(Samples)
    ∫fₓdx = Avg * ∫gₓdx
    Δ = ∫gₓdx * σ / √Ŋ
    return ∫fₓdx, Δ
end

SMCParameters = Dict(:f => x -> ℯ^(- x^2), :X₀ => 0, :X₁ => 2, :Ŋ => 10^7)
IMCParameters = Dict(:f => x -> ℯ^(- x^2), :g => x -> ℯ^-x, :∫gₓdx => 1 - ℯ^-2, :Đₓ => x -> -log(1 - (1 - ℯ^-2) * x), :Ŋ => 10^7)

SMC(;SMCParameters...)
IMC(;IMCParameters...)
