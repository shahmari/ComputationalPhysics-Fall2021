using ProgressBars, Plots, StatsBase, Statistics, LaTeXStrings, JLD

function resultsPart1(N)
    return [rand(0:9) for n ∈ 1:N]
end

function resultsPart2(N)
    TotRandoms = [rand(0:9) for n ∈ 1:N]
    TotRandNumbers = [length(findall(x->x==n,TotRandoms)) for n ∈ 0:9]
    return std(TotRandNumbers)/N
end

function Linear_fit(X, Y)
    A = [hcat(X) reshape(ones(length(X)), length(X), 1)]
    b = reshape(hcat(Y), length(Y), 1)
    line = (A \ b)
    return line
end

N = 10^6
TotRandoms = resultsPart1(N)
histogram(TotRandoms, bins = -0.5:9.5, xticks = 0:9, legend = nothing, framestyle = :box, c = :steelblue)
plot!([-0.5,9.5],[N/10,N/10], linestyle = :dash, c = :black, linewidth = 4)

#=====================##=====================##=====================#


span = ceil.(Int,exp.(2:0.25:10))
Data = zeros(1000,length(span))

for n in ProgressBar(1:1000)
    Data[n,:] = [resultsPart2(N) for N in span]
end

CVList = [mean(Data[:,n]) for n ∈ 1:length(span)]
Line = Linear_fit(log.(span),log.(CVList))
X = log.(span)
Y = X .* Line[1,1] .+ Line[2,1]


plot(X,Y, label = L"Y = %$(round(Line[1],digits=3))\ X + %$(round(Line[2],digits=3))", line = :dash, c = :black)
scatter!(log.(span), log.(CVList),label = L"Data\ Point", framestyle = :box, c = :steelblue)
scatter!(title = L"log-log\ Plot\ CV-N\ (e^2>N>e^{10},\ 10000\ runs)",
    ylabel = L"ln_{(CV_{(N)})}", xlabel = L"ln_{(N)}")
