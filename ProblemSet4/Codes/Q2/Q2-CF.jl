using Plots, StatsBase, Statistics, LaTeXStrings, JLD, ProgressMeter, LsqFit

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


dimlist = [10,15,20,25,30,35,40,45,50,55,60,65,70,75,100,125,150]
# runlist = [10000,5000,2000,1000,500,200,100,100,100]
PList = hcat(0.52:0.001:0.62)
runnum = 100
progress = Progress(length(PList)*runnum*length(dimlist); showspeed=true)
# progress = Progress(length(PList)*sum(runlist); showspeed=true)
for n in 1:length(dimlist)
    dim = dimlist[n]
    AvgData = []
    STDData = []
    TopAvg = []
    for p in PList
        totalruns = []
        for i in 1:runnum
            push!(totalruns, ReturnXi(dim, p))
            next!(progress)
            update!(progress)
        end
        push!(TopAvg,mean(sort(totalruns,rev = true)[1:Int(runnum/5)]))
        push!(AvgData,mean(totalruns))
        push!(STDData,std(totalruns))
    end
    save("../../Data/Q2/Q2-$dim-fc.jld", "TopAvg", TopAvg, "AvgData", AvgData, "STDData", STDData)
end

Data = zeros(length(dimlist),length(PList),3)

for n in 1:length(dimlist)
    print("\r$n")
    dim = dimlist[n]
    data = load("../../Data/Q2/Q2-$dim-fc.jld")
    Data[n,:,1] = data["TopAvg"]
    Data[n,:,2] = data["AvgData"]
    Data[n,:,3] = data["STDData"]
end

PmaxList = []
for n in 1:length(dimlist)
    push!(PmaxList, PList[findfirst(x->x==max(Data[n,:,2]...),Data[n,:,2])])
end

@. modelP∞(x, p) = abs(x - p[1])^(-1.2)
@. modelν(x, p) = abs(x - 0.5927)^(-p[1])

xdata = PmaxList
ydata = dimlist
InitP∞ = [0.59]
Initν = [1.2]

ν = curve_fit(modelν, xdata, ydata, Initν).param[1]
P∞ = curve_fit(modelP∞, xdata, ydata, InitP∞).param[1]


P∞range = hcat(0.52:0.001:0.596)
νrange = abs.(P∞range .- P∞).^(-ν)

scatter(PmaxList,dimlist, label=L"Data\ Points", title=L"Curve\ fitting\ Plot\ (\nu = %$(round(ν,digits= 2)) ,\ Pc_{\infty} = %$(round(P∞,digits= 2)))",legend = 150)
plot!(P∞range,νrange,c= :purple,lw = 2,label = L"L=|Pc_{(L)} - Pc_{\infty}|^{-\nu}")
plot!([0.6,0.6],[0.0,200.0], linestyle = :dash,c = :black,lw = 2,
    label = L"X=Pc_{\infty}= %$(round(P∞,digits= 2))")
savefig("../../Figs/Q2/Q2-CF-res.pdf")
