using Plots, StatsBase, Statistics, LaTeXStrings, JLD

function FindBorder(dim, Network)
    border0 = []
    for i in 1:dim
        for j in 1:dim
            check = false
            if Network[i,j] == 0
                Neighbors = ((i-1,j),(i,j-1),(i+1,j),(i,j+1))
                for nei in Neighbors
                    if dim+1 ∉ nei && 0 ∉ nei && Network[nei...] == 1
                        check = true
                        break
                    end
                end
                if check == true
                    push!(border0,(i,j))
                end
            end
        end
    end
    return border0
end

function FindCorrelationLength(Network,dim)
    FindOnSpots = findall(x->x==1,Network)
    SpotNums = length(FindOnSpots)
    ΣR² = 0.0; iMC = 0.0; jMC = 0.0
    for indx in FindOnSpots
        iMC += indx[1]/SpotNums
        jMC += indx[2]/SpotNums
    end
    for indx in FindOnSpots
        ΣR² += ((indx[1]-iMC)^2 + (indx[2]-jMC)^2)/SpotNums
    end
    return sqrt(ΣR²)
end

function ClusterGrowth(dim,P)
    Network = zeros(Int,dim,dim)
    # Network[rand(1:dim^2)] = 1
    Network[round(Int, dim/2), round(Int, dim/2)] = 1
    while 0 in Network
        Borders = FindBorder(dim,Network)
        if length(Borders) == 0
            break
        end
        for spot in Borders
            if P > rand()
                Network[spot...] = 1
            else
                Network[spot...] = -1
            end
        end
    end
    return Network
end
