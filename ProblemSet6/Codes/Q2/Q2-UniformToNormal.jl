using ProgressBars, Plots, StatsBase, LaTeXStrings, JLD, Distributions, StatsPlots

Plots.default(titlefontsize=12, tickfontsize=10, labelfontsize=16,
    fontfamily="Computer Modern")

function UniformRand(;N,a,c,m)
    seed = rand(1:1000)

    series = [seed]
    for n ∈ 1:N+1
        randvalue = (a*series[end] + c) % m
        push!(series,randvalue)
    end
    return series[3:end]
end

# function reshuffle(parameters)
#     return round.(Int,baserand(;parameters...) * (parameters[:N]-1) / parameters[:m] .+ 1)
# end

function FloatUniformRand(N)
    parameters = Dict(:N => N, :a => 214013, :c => 2531011, :m => 2^32)
    return UniformRand(;parameters...) / parameters[:m]
end

function GaussRandDist(σ, μ, N)
    x₁ = FloatUniformRand(round(Int,N/2))
    x₂ = FloatUniformRand(round(Int,N/2))

    ρ = σ * .√(2 * log.(1 ./ (1 .- x₁)))
    θ = 2π * x₂

    y₁ = ρ .* sin.(θ) .+ μ
    y₂ = ρ .* cos.(θ) .+ μ

    return [y₁; y₂]
end

Data = [GaussRandDist(1,3,1000000),
            GaussRandDist(2,6,1000000),
                GaussRandDist(3,9,1000000),
                    GaussRandDist(4,12,1000000),]

cd(dirname(@__FILE__))

save("../../Data/Q2/Q2-samples.jld","Samples", Data)

colors = [:steelblue, :purple, :green, :brown]
histogram(Data[1],  bins = 25,
    label = L"Generated\ Normal:\ \mu =3 , \sigma =1",
    alpha = 0.5, c = colors[1], normalize = true)
plot!(Normal(3,1), c = :black, label = L"Normal Curve:\ \mathcal{N}(\mu =3 , \sigma =1)")

histogram!(Data[2], bins = 50,
    label = L"Generated\ Normal:\ \mu =6 , \sigma =2",
    alpha = 0.5, c = colors[2], normalize = true)
plot!(Normal(6,2), c = :black, label = L"Normal Curve:\ \mathcal{N}(\mu =6 , \sigma =2)")

histogram!(Data[3], bins = 75,
    label = L"Generated\ Normal:\ \mu =9 , \sigma =3",
    alpha = 0.5, c = colors[3], normalize = true)
plot!(Normal(9,3), c = :black, label = L"Normal Curve:\ \mathcal{N}(\mu =9 , \sigma =3)")

histogram!(Data[4], bins = 100,
    label = L"Generated\ Normal:\ \mu =12 , \sigma =4",
    alpha = 0.5, c = colors[4], normalize = true)
plot!(Normal(12,4), c = :black, label = L"Normal Curve:\ \mathcal{N}(\mu =12 , \sigma =4)")
plot!(title = L"Normalized\ Distribution\ of\ Samples\ (10^6\ sample\ for\ each)",frame = :box)

savefig("../../Figs/Q2/Q2-Hist.pdf")
