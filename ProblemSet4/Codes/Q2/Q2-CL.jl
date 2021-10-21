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
    return Network, S, L
end

function JoinColors(Network, L)
    dim = size(Network)[1]
    for i in 1:dim
        for j in 1:dim
            if Network[i,j]!= -1
                Network[i,j] = RecursionMapping(Network[i,j],L)
            end
        end
    end
    return Network
end

function FindCorrelationLength(Network, S)
    if length(Set(S)) <= 2
        return 0.0
    end
    if findall(x->x==max(S...),S)[1] ∈ intersect(Network[1,:],Network[end,:])
        BiggestFinite = findall(x->x==sort(S, rev=true)[2],S)[1]
    else
        BiggestFinite = findall(x->x==max(S...),S)[1]
    end
    BFiniteCluster = findall(x->x==BiggestFinite,Network)
    BFiniteClusterSize = length(BFiniteCluster)
    dim = size(Network)[1]; ΣR² = 0.0; iMC = 0.0; jMC = 0.0
    for indx in BFiniteCluster
        iMC += indx[1]/BFiniteClusterSize
        jMC += indx[2]/BFiniteClusterSize
    end
    for indx in BFiniteCluster
        ΣR² += ((indx[1]-iMC)^2 + (indx[2]-jMC)^2)/BFiniteClusterSize
    end
    return sqrt(ΣR²)
end

scatter()
dimlist = [10,20,40,80,160]
runnumlist = [10000,5000,3000,2000,1000]
totalData = []
for i in 1:5
    dim = dimlist[i]
    runnum = runnumlist[i]
    CLAvg = []
    CLSTD = []
    PList = hcat(0:0.02:1)
    for p in PList
        totalruns = []
        for i in 1:runnum
            Network, S, L = HKNetworkDynamic(dim, p)
            Network = JoinColors(Network, L)
            push!(totalruns, FindCorrelationLength(Network, S))
        end
        push!(CLAvg, mean(totalruns))
        push!(CLSTD, std(totalruns))
        print("\r$p")
    end
    push!(totalData, (CLAvg,CLSTD))
    scatter!(PList,CLAvg, yerr= CLSTD,label = L"Network\ %$dim\times%$dim,\ %$runnum\ runs")
end
save("../../Data/Q2/Q2-totdata.jld", "data", totalData)
scatter!(
    title = L"Correlation\ Length",
    xlabel = L"P",
    ylabel = L"\xi_{(P)}")
savefig("../../Figs/Q2/CL.pdf")
