using Plots

function ColoringTheNode(i, j, Network, Num)
    if P > rand() && Network[i, j] == 0
        Network[i, j] = Num
        Num += 1
    end
    return Network, Num
end

function FindNeighbors(i, j, Network, dim)
    Numbereds = 0
    colorList = [Network[i, j]]
    for neighbor in ([0,1], [1,0], [0,-1], [-1,0])
        if j + neighbor[2] in 1:dim && Network[([i, j] + neighbor)...] > 0
            push!(colorList,Network[([i, j] + neighbor)...])
            Numbereds += 1
        end
    end
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

P = 0.01
dim = 100
maxnum = 1000

Network = rand((0,-1),dim,dim)

Network[1,:] = fill(1,dim)
Network[end,:] = fill(maxnum,dim)

Num = 2

for i in 2:dim-1
    for j in 1:dim
        Network, Num = ColoringTheNode(i, j, Network, Num)
        colorList, Numbereds = FindNeighbors(i, j, Network, dim)
        if Numbereds == 0
            continue
        elseif Numbereds == 1
            Network[i, j] = colorList[2]
        else
            Network = JoinColors(Network, dim, colorList)
        end
    end
end
