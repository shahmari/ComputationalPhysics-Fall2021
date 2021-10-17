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
    # return CheckPercolation(dim, L, Network)
end

function JoinColors(Network, dim, L)
    for i in 1:dim
        for j in 1:dim
            if Network[i,j]!= -1
                Network[i,j] = RecursionMapping(Network[i,j],L)
            end
        end
    end
    return Network
end

# function FindFinite(Network, dim, L)
#     fside = []
#     lside = []
#     for i in 1:dim
#         if Network[i,dim] != -1
#             push!(lside, RecursionMapping(Network[i,dim],L))
#         end
#         if Network[i,1] != -1
#             push!(fside, RecursionMapping(Network[i,1],L))
#         end
#     end
#     Finites = []
#     for i in L
#         if RecursionMapping(i,L) ∉ intersect(fside,lside) && RecursionMapping(i,L) ∉ Finites
#             push!(Finites, RecursionMapping(i,L))
#         end
#     end
#     return Finites
# end

# function FindCorrelationLength(Network, dim, L)
#     Finites = FindFinite(Network, L)
#     iMC, jMC = FindMassCenter(Network,dim)
#
#     GyrRid = []
#     for val in Finites
#         RList = []
#         for i in 1:dim
#             for j in 1:dim
#                 if Network[i,j] != -1 && RecursionMapping(Network[i,j], L) == val
#                     push!(RList, (iMC - i)^2 + (jMC - j)^2)
#                 end
#             end
#         end
#         push!(GyrRid, sqrt(mean(RList)))
#     end
#     if length(GyrRid) > 0
#         return mean(GyrRid)
#     else
#         return 0.0
#     end
# end

function FindFinite(Network, L)
    Infinites = intersect(Network[1,:],Network[dim,:])
    Finites = []
    for i in L
        if RecursionMapping(i,L) ∉ Infinites && RecursionMapping(i,L) ∉ Finites
            push!(Finites, RecursionMapping(i,L))
        end
    end
    return Finites
end

function FindMassCenter(Network, dim)
    ilist = 0
    jlist = 0
    itnum = 0
    for i in 1:dim
        for j in 1:dim
            if Network[i,j] != -1
                ilist += i
                jlist += j
                itnum += 1
            end
        end
    end
    return ilist/itnum , jlist/itnum
end

function FindCorrelationLength(Network, dim, L)
    Finites = FindFinite(Network, L)
    iMC, jMC = FindMassCenter(Network,dim)

    GyrRid = []
    for val in Finites
        RList = []
        for i in 1:dim
            for j in 1:dim
                if Network[i,j] != -1 && Network[i,j] == val
                    push!(RList, (iMC - i)^2 + (jMC - j)^2)
                end
            end
        end
        push!(GyrRid, sqrt(mean(RList)))
    end
    if length(GyrRid) > 0
        return mean(GyrRid)
    else
        return 0.0
    end
end

# dim = 10
# CLlist = []
# PList = hcat(0:0.02:1)
# for p in PList
#     Network, S, L = HKNetworkDynamic(dim, p)
#     push!(CLlist, FindCorrelationLength(Network, dim, L))
# end
#
# scatter(PList,CLlist)
