using Plots, StatsBase, Statistics

function ColoringTheNode(i, j, Network, Num)
    check = false
    if P > rand() && Network[i, j] == 0
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

function InitialNetwork(dim, Pnet)
    Network = sample([-1500,0], Weights([1-Pnet, Pnet]),(dim,dim))
    return Network
end

function PercolationCheck(dim, Pnet, P)
    Num = 2
    maxnum = 10000
    Network = InitialNetwork(dim, Pnet)
    Network[1,:] = fill(1,dim)
    Network[end,:] = fill(maxnum,dim)
    TotalNets = []
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
            push!(TotalNets,copy(Network))
        end
    end
    return TotalNets
end


P=0.85
save("Percolation$P.jld", "data", PercolationCheck(100,P,1))

Data = PercolationCheck(100,P,1)
anim = Animation()
for i in 110:100:length(Data)
    plt = heatmap(Data[i], c = cgrad(:copper, 10), legend = false, border = :none, title = "Percolation for P = $P")
    frame(anim, plt)
end
gif(anim, "../../Figs/Q4/Percolation$P.gif")
