using Plots, LaTeXStrings, Statistics

function Deposition(;len, tot_time, rate)
    surf = [0 for i=1:len]
    meanList = [0.0 for i=1:tot_time]
    for n in 1:tot_time
        randsurf = rand(1:len,rate)
        for i in randsurf
            surf[i] += 1
        end
        meanList[n] = mean(surf)
    end
    return meanList
end

function Linear_fit(;len, tot_time, rate)
    A = [hcat(1:tot_time) reshape(ones(tot_time), tot_time, 1)]
    b = reshape(meanAvg, tot_time, 1)
    line = (A \ b)
    x = 0:tot_time
    y = x .* line[1] .+ line[2]
    return x, y, line
end

iternum = 1000
Parameters = Dict(
                :len => 200,
                :tot_time => 20,
                :rate => 5000)
allmean = [[0.0 for i in 1:Parameters[:tot_time]] for j = 1:iternum]
meanAvg = [0.0 for i in 1:Parameters[:tot_time]]
means = [0.0 for i in 1:Parameters[:tot_time]]
for i in 1:iternum
    MeanList = Deposition(;Parameters...)
    allmean[i] = MeanList
    meanAvg += MeanList
    print("\r$i")
end
meanAvg /= iternum
for i in 1:Parameters[:tot_time]
    means[i] = std(log.(hcat(allmean...))[i,:])
end

X, Y, Line = Linear_fit(;Parameters...)

scatter(meanAvg,
    xlabel= L"Time",
    ylabel= L"H_{(t)}",
    title= L"H_{(t)}-Time, Particles \ number = 10^{%$(Int(log10(Parameters[:tot_time]*Parameters[:rate])))}",
    yerror = means,
    label = L"Data\ point")
plot!(X,Y,c= :black,label = L"y = %$(round(Line[1],digits= 2))x + %$(round(Line[2],digits= 2))")
savefig("C:\\Users\\Yaghoub\\Documents\\GitHub\\ComputationalPhysics-Fall2021\\ProblemSet2\\Figs\\Q2\\RBDH(t).png")
