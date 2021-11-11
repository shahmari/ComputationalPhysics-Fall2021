using ProgressBars, Plots, StatsBase, LaTeXStrings, JLD, Distributions, StatsPlots

Plots.default(titlefontsize=12, tickfontsize=10, labelfontsize=16,
    fontfamily="Computer Modern")

SampleNumber = 100000
NList = [5, 10, 100, 1000]
Data = [[sum(rand(N)) for i in ProgressBar(1:SampleNumber)] for N ∈ NList]

cd(dirname(@__FILE__))
save("../../Data/Q1/Q1-samples.jld","Samples", Data)

plots = []
colors = [:steelblue, :purple, :green, :brown]
for i ∈ 1:4
    Plt = begin
    histogram(Data[i],nbins = 60, normalize = true, c = colors[i], label = L"Histogram\ of\ Samples")
    plot!(Normal(mean(Data[i]),std(Data[i])), c = :black, lineWidth = 10,
        label = L"\mathcal{Normal}(\mu = %$(round(mean(Data[i]),digits = 2)), \sigma = %$(round(std(Data[i]),digits = 2)))")
    plot!(frame = :box, title = "N = $(NList[i]), $SampleNumber Samples",
        xlabel = L"X_N", ylabel = L"Counting\ ratio")
    end
    push!(plots,Plt)
end


plot(plots..., size = (1000,1000), plot_title = "Normalized Distribution of Samples")
savefig("../../Figs/Q1/Q1-Hist.pdf")
