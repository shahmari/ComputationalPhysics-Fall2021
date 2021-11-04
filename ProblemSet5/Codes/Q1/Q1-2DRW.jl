using ProgressBars, Plots, StatsBase, Statistics, LaTeXStrings, JLD

@inline function RandomWalker2D(TotalSteps, R₀)
    choices = ((1,0),(0,1),(-1,0),(0,-1))
    TotalWalks = rand(choices,TotalSteps)
    return reduce(.+,TotalWalks) .+ R₀
end


function RGyration(ensemble,RWNumbers)
    TotalR = collect(Tuple(ensemble[i,:]) for i ∈ 1:RWNumbers)
    MassCenter = (mean(TotalR[1]),mean(TotalR[2]))
    return √(sum(pos -> reduce(.+, (pos .- MassCenter).^2), TotalR)/RWNumbers)
end
#=
function PlotHeatmap(Data, RWNumbers) #useless piece of shit
    Tuples = collect(Tuple(Data[i,:]) for i ∈ 1:RWNumbers)
    maxindx = max(getindex.(Tuples, 1)...,getindex.(Tuples, 2)...)
    minindx = min(getindex.(Tuples, 1)...,getindex.(Tuples, 2)...)
    delta = maxindx - minindx + 1
    Surf = zeros(Int,delta,delta)

    for tup in Tuples
        Surf[(tup .- (minindx - 1))...] += 1
    end
    return heatmap(Surf)
end
=#

TStepsList = [10,100,1000]
RWNumbers = 100000
R₀ = [0,0]
Data = []
for TotalSteps in TStepsList
    TotalWalkers = zeros(Int,RWNumbers,2)
    for i ∈ ProgressBar(1:RWNumbers)
        TotalWalkers[i,:] = RandomWalker2D(TotalSteps,R₀)
    end
    push!(Data, copy(TotalWalkers))
end
save("../../Data/Q1/Q1-2DHist.jld", "Data", Data)

#Plot section

scatter(RgList.^2)

# PlotHeatmap(Data, RWNumbers)
dimnum = 3
P1 = histogram2d(collect(Tuple(Data[dimnum][i,:]) for i ∈ 1:RWNumbers),
    nbin = 200,
    xaxis = nothing,
    yaxis = nothing,
    background=:black,
    framestyle = :box
    )
P2 = histogram(Data[dimnum][:,2], orientation = :horizontal, xticks = [0,3000],
    ylabel = L"Y", flip = true,framestyle = :box, c = :purple)
P3 = histogram(Data[dimnum][:,1], flip = true, yticks = [0,3000],
    xlabel = L"X", framestyle = :box, c = :purple)
P4 = plot(legend=false,grid=false,
    foreground_color_subplot=:black,background=:black,framestyle = :box)

P21 = plot(P2,P1,layout = @layout[a{0.1w} grid(1,1, heights = [0.95])], legend = nothing)

P22 = plot(P4,P3,P4,layout = @layout[a{0.08w} grid(1,1) a{0.02w}], legend = nothing)

plot(P21,P22,layout = grid(2,1, heights = [0.85,0.15]),
    plot_title = L"Distribution\ of\ %$RWNumbers\ RW\ in\ %$(TStepsList[dimnum])\ steps")

savefig("../../Figs/Q1/Q1-2DHist$dimnum.pdf")
