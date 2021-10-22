using Plots, StatsBase, Statistics, LaTeXStrings, JLD

dim = 100

Network = zeros(Int,dim,dim)
Network[rand(1:dim^2)] = 1

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

FindBorder(dim,Network)
