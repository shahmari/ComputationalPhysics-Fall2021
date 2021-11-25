using Plots, Distributions, StatsPlots, LaTeXStrings, StatsBase



function Metropolis(P::Function, x₀::Real = 0.0, Δ::Real = 2.9, Steps::Integer = 1000000)
    Aₙ = 0
    XList = Vector{Float64}(undef, Steps)
    x = x₀
    for n = 1:Steps
        XList[n] = x
        y = x + rand(Uniform(-Δ, Δ))
        if rand() < P(y) / P(x)
            x = y
            Aₙ += 1
        end
    end
    return XList, Aₙ / Steps
end

Data = []

histogram(Metropolis(x -> ℯ^(-x^2 / 2))[1], bins = 100, normalize = true)
plot!(Normal(0, 1))

Metropolis(x -> ℯ^(-x^2 / 2))[2]

Data = Metropolis(x -> ℯ^(-x^2 / 2))[1]

autocor(Metropolis(x -> ℯ^(-x^2 / 2))[1])