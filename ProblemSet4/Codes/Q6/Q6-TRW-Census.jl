using ProgressMeter, Plots, StatsBase, Statistics, LaTeXStrings, JLD

function TRW_CD(L,P,X0)
    surf = zeros(L+2); surf[X0] = 1.0
    TimePassed = 0
    while surf[1] + surf[end] < 0.99
        Altsurf = copy(surf)
        for i ∈ 2:L+1
            Altsurf[i+1] += surf[i]*(1-P)
            Altsurf[i-1] += surf[i]*P
            Altsurf[i] -= surf[i]
        end
        surf = Altsurf
        TimePassed += 1
    end
    return TimePassed
end

L = 40
PList = hcat(0:0.01:1)
X0List = hcat(1:L+2)

Data = zeros(length(X0List),length(PList))
progress = Progress(length(X0List)*length(PList)*runnum)

for i in 1:length(X0List)
    X0 = X0List[i]
    for j in 1:length(PList)
        P = PList[j]
        Data[i,j] = mean(TRW_CD(L,P,X0))
    end
end

save("../../Data/Q6/Q6-LT.jld", "Data", Data)
Data = load("../../Data/Q6/Q6-LT.jld")["Data"]

P1 = heatmap(0:0.01:1, 1:L+2,Data, xaxis=nothing,legend = nothing,title = L"Life\ time")
P2 = contour(0:0.01:1, 1:L+2,Data, fill = true ,legend = nothing,xlabel = L"P")
P3 = plot(P1,P2,layout = (2,1),ylabel = L"X_{init}",link = :y)
CBar = heatmap((0:1650).*ones(1651,1),xlabel ="",legend=:none, xticks=:none, yticks=(1:165:1651, string.(0:165:1651)))

plot([P3]...,CBar,layout = @layout[grid(1,1) a{0.04w}], size = (600,500))
savefig("../../Figs/Q6/Q6-LT.pdf")

X0List[21]


X0Plot = []
for X0 in [11,21,31]
    PLT = begin
        plot(PList,Data[X0,:,1], c = :black,linewidth = 3)
        plot!(PList,Data[X0,:,1], linewidth = 2, c = :purple, fill = (0, 0.5))
        plot!(title = L"X_{init} = %$(X0List[X0])",
            xlabel = L"P", ylabel = L"Avg\ Life\ Time", legend = :none)
    end
    push!(X0Plot, PLT)
end

X0P = plot(X0Plot...,layout = (1,3))

PPlot = []
for P in [26,51,76]
    PLT = begin
        plot(X0List,Data[:,P], c = :black, linewidth = 3)
        plot!(X0List,Data[:,P], linewidth = 2, c = :steelblue,  fill = (0, 0.5))
        plot!(title = L"P = %$(PList[P])",
            xlabel = L"X0", ylabel = L"Avg\ Life\ Time", legend = :none)
    end
    push!(PPlot, PLT)
end
PPlot[1]
PP = plot(PPlot...,layout = (1,3))

Pmain2 = plot(PP,X0P,layout = (2,1),size = (850,690), plot_title = L"Life\ time\ of\ a\ Random\ Walker\ (%$runnum \ runs)")
savefig(Pmain2,"../../Figs/Q6/Q6-scat.pdf")
