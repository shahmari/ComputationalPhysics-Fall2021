using Plots, StatsBase, Statistics, LaTeXStrings, JLD, ProgressMeter

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
    if findfirst(x->x==max(S...),S) ∈ intersect(Network[1,:],Network[end,:])
        BiggestFinite = findfirst(x->x==sort(S, rev=true)[2],S)
    else
        BiggestFinite = findfirst(x->x==max(S...),S)
    end
    BFiniteCluster = Tuple.(findall(x->x==BiggestFinite,Network))
    BFiniteClusterSize = length(BFiniteCluster)
    MassCenter = reduce(.+, BFiniteCluster) ./ BFiniteClusterSize
    return √(sum(pos -> reduce(.+, (pos .- MassCenter).^2),Tuple.(BFiniteCluster))/BFiniteClusterSize)
end

function ReturnXi(dim, p)
    Network, S, L = HKNetworkDynamic(dim, p)
    Network = JoinColors(Network, L)
    return  FindCorrelationLength(Network, S)
end


dimlist = [10,20,30,40,50,75,100,125,150]
PList = hcat(0.5:0.005:0.65)
progress = Progress(length(PList)*length(dimlist))
for dim in dimlist
    runnum = 2500000/ dim^2
    AvgData = []
    for p in PList
        totalruns = []
        for i in 1:runnum
            push!(totalruns, ReturnXi(dim, p))
            print("\r$i    ")
        end
        push!(AvgData,mean)
        next!(progress)
    end
    print("\r$n    ")
    save("../../Data/Q2/Q2-$dim-fc.jld", "data", totalruns)
end

progress = Progress(300)
for i in 1:3
    for j in 1:100
        sleep(0.075)
        next!(progress)
    end
    sleep(0.5)
end
