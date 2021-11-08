using ProgressBars, Plots, StatsBase, Statistics, LaTeXStrings, JLD

choices = ((1,0),(0,1),(-1,0),(0,-1))

function RecurseSAW(R, Network, N)
    if N == 0
        return 1
    end
    count_ = 0
    Network[R...] = true
    for step in choices
        Rₙ = R .+ step
        if !Network[Rₙ...]
            count_ += RecurseSAW(Rₙ, Network, N-1)
        end
    end
    Network[R...] = false
    return count_
end

function maincal(N)
    dim = 2N+1
    Network = falses(dim,dim)
    R = [N+1,N+1]
    return RecurseSAW(R, Network, N)
end

maincal(15)
