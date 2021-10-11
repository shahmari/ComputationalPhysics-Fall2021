using Plots, StatsBase

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

P = 1
Pnet = 0.8
dim = 100
maxnum = 100000

Network = sample([-1,0], Weights([1-Pnet, Pnet]),(dim,dim))

Network[1,:] = fill(1,dim)
Network[end,:] = fill(maxnum,dim)

Num = 2

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

Network
