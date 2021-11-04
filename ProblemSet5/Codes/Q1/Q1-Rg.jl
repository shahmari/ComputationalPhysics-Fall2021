using ProgressBars, Plots, StatsBase, Statistics, LaTeXStrings, JLD

function RandomWalker2D(TotalSteps)
    choices = ((1,0),(0,1),(-1,0),(0,-1))
    return reduce(.+,rand(choices,TotalSteps))
end

function RGyration(ensemble)
    MassCenter = (mean(getindex.(ensemble,1)),mean(getindex.(ensemble,2)))
    return √(sum(pos -> reduce(.+, (pos .- MassCenter).^2), ensemble)/length(ensemble))
end

TotalSteps = 40
RWNumbers = 10000000

Data = [[RandomWalker2D(TS) for i ∈ 1:RWNumbers] for TS ∈ ProgressBar(1:TotalSteps)]

RgList = [RGyration(Data[i]) for i ∈ 1:TotalSteps]
save("../../Data/Q1/Q1-Rg.jld", "RgList", RgList)


plot(1:TotalSteps,sqrt.(1:TotalSteps), label = L"Y = \sqrt{X}", c = :black, linestyle = :dash)
scatter!(RgList, label = L"Data\ Point", framestyle = :box, c = :steelblue, legend = 140)
plot!(xlabel = L"Time", ylabel = L"{R_g}_{(t)}",
    title = L"Radius\ of\ Gyration\ of\ RW\ per\ Time\ for\ 10^6 \ Particles")
savefig("../../Figs/Q1/Q1-Rg(t).pdf")
