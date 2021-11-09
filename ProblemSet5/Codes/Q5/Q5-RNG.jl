using ProgressBars, Plots, StatsBase, Statistics, LaTeXStrings, JLD

function findseries(N)
    TotalRandoms = rand(0:9,12*N)
    selected = findall(x-> TotalRandoms[x+1]==4,1:12N-1)
    return TotalRandoms[selected][1:N]
end

function findseries2(N)
    TotalRandoms = []
    while length(TotalRandoms) != N
        randomnums = rand(0:9,2)
        if randomnums[2] == 4
            push!(TotalRandoms,randomnums[1])
        end
    end
    return TotalRandoms
end

function CVCal(N)
    TotRandoms = findseries2(N)
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
TotRandoms = findseries(N)
save("../../Data/Q5/Q5-hist.jld", "Data", TotRandoms)

plot([-0.5,9.5],[N/10,N/10], linestyle = :dash, label = L"Y = \frac{N}{10}", c = :black, linewidth = 4, ylims = (-5000,1.35*10^5))
histogram!(TotRandoms, bins = -0.5:9.5, xticks = 0:9, label = L"Histogram", framestyle = :box, c = :steelblue)
plot!(xlabel = L"Number", ylabel = L"Count",title = L"%$N\ Random\ Numbers\ in\ [0, 9]")
savefig("../../Figs/Q5/Q5-Hist.pdf")

#=====================##=====================##=====================#


span_ = ceil.(Int,exp.(2:0.25:10))
Data = zeros(1000,length(span_))

for n in ProgressBar(1:1000)
    Data[n,:] = [CVCal(N) for N in span_]
end

CVList = [mean(Data[:,n]) for n ∈ 1:length(span_)]
Line = Linear_fit(log.(span_),log.(CVList))
X = log.(span_)
Y = X .* Line[1,1] .+ Line[2,1]

save("../../Data/Q5/Q5-CV.jld", "Data", CVList)

plot(X,Y, label = L"Y = %$(round(Line[1],digits=3))\ X + %$(round(Line[2],digits=3))", line = :dash, c = :black)
scatter!(log.(span_), log.(CVList),label = L"Data\ Point", framestyle = :box, c = :steelblue)
scatter!(title = L"log-log\ Plot\ CV-N\ (e^2>N>e^{10},\ 10000\ runs)",
    ylabel = L"ln_{(CV_{(N)})}", xlabel = L"ln_{(N)}")
savefig("../../Figs/Q5/Q5-Scat.pdf")
