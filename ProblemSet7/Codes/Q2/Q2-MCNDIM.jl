using Plots, Statistics, LaTeXStrings, JLD

function SMC(; f::Function, R₀::Real = 0.0, R₁::Real, θ₀::Real = 0.0, θ₁::Real = π, Ŋ::Integer = 10^5)
    ΔR = R₁ - R₀
    Δθ = θ₁ - θ₀
    Đist = [
        muladd.(rand(Ŋ), ΔR, R₀),
        muladd.(rand(Ŋ), Δθ, θ₀)]
    fₓ = [f(getindex.(Đist, n)...) for n ∈ 1:Ŋ]
    Avg = mean(fₓ)
    ∫fₓdx = Avg * ΔR * Δθ
    Δ = std(fₓ) / √Ŋ
    return ∫fₓdx, Δ
end

SMCParameters = Dict(:f => (r, θ) -> r^3 * sin(θ) * cos(θ) * (3 + r * cos(θ)), :R₁ => 1)

SMC(; SMCParameters...)
