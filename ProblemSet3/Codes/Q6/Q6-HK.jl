using Plots, StatsBase, Statistics, LaTeXStrings, JLD

function InitialNetwork(dim, P)
    Network = sample([-1,0], Weights([1-P, P]),(dim,dim))
    return Network
end

function FindNeighbors(i, j, Network, dim)
    Numbereds = 0
    colorList = []
    for neighbor in ([0,-1], [-1,0])
        if j + neighbor[2] >= 1 && i + neighbor[1] >= 1 && Network[([i, j] + neighbor)...] > 0
            push!(colorList,Network[([i, j] + neighbor)...])
            Numbereds += 1
        end
    end
    return colorList, Numbereds #colorList contain the node and its neighbors
end

function RecursionMapping(x, L)
    while L[x] != L[L[x]]
        x = L[x]
    end
    return L[x]
end

@inline function CheckPercolation(dim, L, Network)
    fside = []
    lside = []
    for i in 1:dim
        if Network[i,dim] != -1
            push!(lside, RecursionMapping(Network[i,dim],L))
        end
        if Network[i,1] != -1
            push!(fside, RecursionMapping(Network[i,1],L))
        end
    end
    if length(intersect(fside,lside)) > 0
        return 1
    else
        return 0
    end
end

@inline function HKNetworkDynamic(dim, P)
    Network = InitialNetwork(dim, P)
    S = []
    L = []
    k = 1

    for j in 1:dim
        for i in 1:dim
            if Network[i,j] != 0
                continue
            end
            NColor, NNum = FindNeighbors(i,j,Network,dim)
            if NNum == 0
                Network[i,j] = k
                push!(L, k)
                push!(S, 1)
                k += 1
            elseif NNum == 1
                Network[i,j] = RecursionMapping(NColor[1],L)
                S[RecursionMapping(NColor[1],L)] += 1
            elseif NNum == 2 && RecursionMapping(NColor[1],L) == RecursionMapping(NColor[2],L)
                Network[i,j] = RecursionMapping(NColor[1],L)
                S[RecursionMapping(NColor[1],L)] += 1
            else
                S[RecursionMapping(NColor[1],L)] += 1 + S[RecursionMapping(NColor[2],L)]
                S[RecursionMapping(NColor[2],L)] = 0
                Network[i,j] = RecursionMapping(NColor[1],L)
                L[NColor[2]] = RecursionMapping(NColor[1],L)
            end
        end
    end
    # return Network
    return CheckPercolation(dim, L, Network)
end

@time HKNetworkDynamic(100, 0.6)


Plist = hcat(0.55:0.002:0.65)
Avglist = []
STDlist = []
TotData = []
# PercolationCheck(40,0.55,1)
for P in Plist
    print("\r$P")
    AllRuns = []
    for i in 1:1000
        push!(AllRuns, HKNetworkDynamic(200, P))
    end
    push!(TotData, AllRuns)
    push!(STDlist, std(AllRuns))
    push!(Avglist, mean(AllRuns))
end

save("200x200-HK.jld", "data", TotData)

scatter(Plist,
    Avglist,
    yerr = STDlist,
    legend = nothing,
    xlabel = L"P",
    ylabel = L"Avg\ Q",
    title = L"Percolation \ for \ 200\times200 \ Network\ (1000 \ runs)")
savefig("../../Figs/Q5/200x200.pdf")
