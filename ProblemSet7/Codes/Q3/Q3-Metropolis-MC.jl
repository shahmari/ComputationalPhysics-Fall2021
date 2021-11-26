module Metropolis

export Metropolis

using Distributions

function Metropolis(; P::Function, x₀::Real = 0.0, Δ::Real = 1.0, Steps::Integer = 1000000)
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

end