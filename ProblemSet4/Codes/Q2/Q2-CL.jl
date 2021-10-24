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

dimlist = [10,20,40,80,160]
runnumlist = [10000,5000,3000,2000,1000]
PList = hcat(0:0.02:1)
totalData = load("../../Data/Q2/Q2-totdata.jld")["data"]

p1 = begin
    scatter(dpi = 200)
    for i in 1:5
        dim = dimlist[i]
        runnum = runnumlist[i]
        data = totalData[i]
        plot!(PList,data[1],linestyle = :dash, label = nothing, c = :black)
        scatter!(PList,data[1],label = L"Network\ %$dim\times%$dim,\ %$runnum\ runs",markersize = 3.5)
    end
    scatter!(
        title = L"Correlation\ Length",
        xlabel = L"P",
        ylabel = L"\xi_{(P)}")
end

p2 = begin
    scatter(dpi = 200)
    for i in 1:5
        dim = dimlist[i]
        runnum = runnumlist[i]
        data = totalData[i]
        scatter!(PList,data[1], ribbon= data[2],markersize = 1)
        plot!(PList,data[1],legend = nothing, c = :black,linestyle = :dot)
    end
    scatter!(
        title = L"Correlation\ Length",
        xlabel = L"P",
        ylabel = L"\xi_{(P)}")
end

p3 = begin
    plot(dpi = 200)
    for i in 1:5
        dim = dimlist[i]
        runnum = runnumlist[i]
        data = totalData[i]
        scatter!(PList,data[1],label = nothing,markersize = 3, c = :black, alpha = 0.5)
        plot!(PList,data[1], label = L"Network\ %$dim\times%$dim,\ %$runnum\ runs")
    end
    scatter!(
        title = L"Correlation\ Length",
        xlabel = L"P",
        ylabel = L"\xi_{(P)}")
end

totalData = load("../../Data/Q2/Q2-totdata-zoomed.jld")["data"]
dimlist = [80,160]
runnumlist = [2000,1000]
PList = hcat(0.5:0.005:0.7)

p4 = begin
    plot()
    col = [:purple,:steelblue]
    for i in 1:2
        dim = dimlist[i]
        runnum = runnumlist[i]
        data = totalData[i]
        plot!(PList,data[1],c = :black, label = nothing)
        scatter!(PList,data[1],ribbon =  data[2],yerr = data[2],label = L"Network\ %$dim\times%$dim,\ %$runnum\ runs", c = col[i])
    end
    scatter!(
        title = L"Correlation\ Length",
        xlabel = L"P",
        ylabel = L"\xi_{(P)}")
end
savefig("../../Figs/Q2/CL3.pdf")

plot(plot(p2,p1, layout= (2,1)),p4,size = (1200,750))
savefig("../../Figs/Q2/CL.pdf")

# for n in 1:5
#     dim = dimlist[n]
#     runnum = runnumlist[n]
#     CLAvg = []
#     CLSTD = []
#     for p in PList
#         totalruns = []
#         for i in 1:runnum
#             Network, S, L = HKNetworkDynamic(dim, p)
#             Network = JoinColors(Network, L)
#             push!(totalruns, FindCorrelationLength(Network, S))
#         end
#         push!(CLAvg, mean(totalruns))
#         push!(CLSTD, std(totalruns))
#         print("\r$p")
#     end
#     push!(totalData, (CLAvg,CLSTD))
# end
# save("../../Data/Q2/Q2-totdata.jld", "data", totalData)

# for n in 1:2
#     dim = dimlist[n]
#     runnum = runnumlist[n]
#     CLAvg = []
#     CLSTD = []
#     for p in PList
#         totalruns = []
#         for i in 1:runnum
#             Network, S, L = HKNetworkDynamic(dim, p)
#             Network = JoinColors(Network, L)
#             push!(totalruns, FindCorrelationLength(Network, S))
#         end
#         push!(CLAvg, mean(totalruns))
#         push!(CLSTD, std(totalruns))
#         print("\r$p")
#     end
#     push!(totalData, (CLAvg,CLSTD))
# end
# save("../../Data/Q2/Q2-totdata-zoomed.jld", "data", totalData)
