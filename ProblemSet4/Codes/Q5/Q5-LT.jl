using ProgressMeter, Plots, StatsBase, Statistics, LaTeXStrings, JLD

function RandomWalker(P,L,X0)
    X = X0
    Tₛ= 0
    while abs(X) != round(Int,L)
        X += sample([1,-1], Weights([P,1-P]))
        Tₛ+= 1
    end
    return Tₛ
end

PList = hcat(0:0.025:1)
X0List = hcat(-20:20)
runnum = 1000

Data = zeros(length(X0List),length(PList),2)
progress = Progress(length(X0List)*length(PList)*runnum)

for i in 1:length(X0List)
    X0 = X0List[i]
    for j in 1:length(PList)
        P = PList[j]
        stepruns = []
        for run in 1:runnum
            push!(stepruns,RandomWalker(P,20,X0))
            next!(progress)
            update!(progress)
        end
        Data[i,j,1] = mean(stepruns)
        Data[i,j,2] = std(stepruns)
    end
end

save("../../Data/Q5/Q5-LT.jld", "Data", Data)

P1 = heatmap(0:0.025:1,-20:20,Data[:,:,1], xaxis=nothing,legend = nothing, yaxis=nothing,title = L"Average\ Life\ time")
P2 = contour(0:0.025:1,-20:20,Data[:,:,1], xaxis=nothing, fill = true ,legend = nothing,ylabel = L"X_{init}",title = L"Average\ Life\ time")
P3 = heatmap(0:0.025:1,-20:20,Data[:,:,2], yaxis=nothing, legend = nothing, xlabel = L"P",title = L"STD\ of\ Life\ time")
P4 = contour(0:0.025:1,-20:20,Data[:,:,2], fill = true ,legend = nothing, xlabel = L"P",ylabel = L"X_{init}",title = L"STD\ of\ Life\ time")

Pmain1 = begin
    P_tot = [P2 P4;P1 P3]
    CBar = plot(heatmap((0:400).*ones(401,1),xlabel ="",
        legend=:none, xticks=:none, yticks=(1:40:401, string.(0:40:400))),
            heatmap((0:300).*ones(301,1),xlabel = L"range",
                legend=:none, xticks=:none, yticks=(1:30:301, string.(0:30:300))),
                    layout = (2,1))
    plot(P_tot...,CBar,
        layout = @layout[grid(2,2) a{0.035w}],
        plot_title = L"Life\ time\ of\ a\ Random\ Walker",size = (800,600))
end
savefig(Pmain1,"../../Figs/Q5/Q5-cont.pdf")



X0Plot = []
for X0 in [11,21,31]
    PLT = begin
        plot(PList,Data[X0,:,1], linestyle = :dash, c = :black)
        scatter!(PList,Data[X0,:,1], yerr = Data[X0,:,1], c = :black, markersize = 5)
        scatter!(PList,Data[X0,:,1], ribbon = Data[X0,:,1], markersize = 3, c = :purple)
        plot!(title = L"X_{init} = %$(X0List[X0])",
            xlabel = L"P", ylabel = L"Avg\ Life\ Time", legend = :none)
    end
    push!(X0Plot, PLT)
end

X0P = plot(X0Plot...,layout = (1,3))

PPlot = []
for P in [11,21,31]
    PLT = begin
        plot(X0List,Data[P,:,1], linestyle = :dash, c = :black)
        scatter!(X0List,Data[P,:,1], yerr = Data[P,:,1], c = :black, markersize = 5)
        scatter!(X0List,Data[P,:,1], ribbon = Data[P,:,1], markersize = 3, c = :steelblue)
        plot!(title = L"P = %$(PList[P])",
            xlabel = L"X0", ylabel = L"Avg\ Life\ Time", legend = :none)
    end
    push!(PPlot, PLT)
end

PP = plot(PPlot...,layout = (1,3))

Pmain2 = plot(PP,X0P,layout = (2,1),size = (850,690), plot_title = L"Life\ time\ of\ a\ Random\ Walker\ ")
savefig(Pmain2,"../../Figs/Q5/Q5-scat.pdf")
