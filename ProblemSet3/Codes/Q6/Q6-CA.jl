using Plots, StatsBase, Statistics, LaTeXStrings, JLD

function ColoringTheNode(i, j, Network, Num)
    check = false
    if Network[i, j] == 0
        Network[i, j] = Num
        check = true
        Num += 1
    end
    return Network, Num, check
end

function FindNeighbors(i, j, Network, dim)
    Numbereds = 0
    colorList = []
    for neighbor in ([0,1], [1,0], [0,-1], [-1,0])
        if j + neighbor[2] >= 1 && j + neighbor[2] <= dim && Network[([i, j] + neighbor)...] > 0
            push!(colorList,Network[([i, j] + neighbor)...])
            Numbereds += 1
        end
    end
    push!(colorList, Network[i, j])
    return colorList, Numbereds #colorList contain the node and its neighbors
end

function JoinColors(Network, dim, colorList)
    mincolor = min(colorList...)
    for i in 1:dim
        for j in 1:dim
            if Network[i, j] in colorList
                Network[i, j] = mincolor
            end
        end
    end
    return Network
end

function InitialNetwork(dim, P)
    Network = sample([-1,0], Weights([1-P, P]),(dim,dim))
    return Network
end

@inline function PercolationCheck(dim, P)
    Num = 2
    maxnum = 100000
    Network = InitialNetwork(dim, P)
    Network[1,:] = ones(dim)
    Network[end,:] = fill(maxnum,dim)
    for i in 2:dim-1
        for j in 1:dim
            Network, Num, check = ColoringTheNode(i, j, Network, Num)
            if check == true
                colorList, Numbereds = FindNeighbors(i, j, Network, dim)
                if Numbereds == 0
                    continue
                elseif Numbereds == 1
                    Network[i, j] = colorList[1]
                else
                    Network = JoinColors(Network, dim, colorList)
                end
            end
        end
    end
    if 1 in Network[end,:]
        return 1
    else
        return 0
    end
end

a = time()
PercolationCheck(100,0.6)
time()-a
Plist = hcat(0:0.02:1)
Avglist = []
STDlist = []
TotData = []
# PercolationCheck(40,0.55,1)
for P in Plist
    print("\r$P")
    AllRuns = []
    for i in 1:5000
        push!(AllRuns, PercolationCheck(10,P))
    end
    push!(TotData, AllRuns)
    push!(STDlist, std(AllRuns))
    push!(Avglist, mean(AllRuns))
end

save("10x10.jld", "data", TotData)

scatter(Plist, Avglist, yerr = STDlist, legend = nothing, xlabel = L"P", ylabel = L"Avg\ Q", title = L"Percolation \ for \ 10\times10 \ Network\ (5000 \ runs)")
savefig("../../Figs/Q5/10x10.pdf")
