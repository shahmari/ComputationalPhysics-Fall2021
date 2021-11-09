using ProgressBars, Plots, LaTeXStrings, JLD

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

TotList = [maincal(n) for n ∈ ProgressBar(0:17)]
save("../../Data/Q3/Q3-SAW.jld", "Data", TotList)

P1 = begin
    plot(0:17,log.(TotList), c = :black, label = nothing)
    scatter!(0:17,log.(TotList), c = :steelblue,framestyle = :box,
        ylabel = L"ln\ (Number\ of\ Paths)", xlabel = L"N",
        label = L"Self-Avoiding\ Walker", legend = 140)
        plot!(0:17, log.(4 .^(0:17)), c = :black, label = nothing)
        scatter!(0:17, log.(4 .^(0:17)), c = :purple, label = L"Random\ Walker")
end


P2 = begin
    plot(0:17,TotList, c = :black, label = nothing)
    scatter!(0:17,TotList, c = :steelblue,framestyle = :box,
        ylabel = L"Number\ of\ Paths", xlabel = L"N",
        label = L"Self-Avoiding\ Walker", legend = 140)
        plot!(0:13,4 .^(0:13), c = :black, label = nothing)
        scatter!(0:13, 4 .^(0:13), c = :purple, label = L"Random\ Walker")
end

savefig(P1,"../../Figs/Q3/Q3-1.pdf")
savefig(P2,"../../Figs/Q3/Q3-2.pdf")

P3 = begin
    plot(0:17, TotList./ 4 .^(0:17), c = :black, label = nothing)
    scatter!(0:17, TotList./ 4 .^(0:17), c = :steelblue,framestyle = :box,
        ylabel = L"Ratio\ of\ Number\ of\ Paths", xlabel = L"N", legend = nothing,
            title = L"Ratio\ of\ Number\ of\ Paths\ of\ SAW\ and\ RW")
end

savefig(P3,"../../Figs/Q3/Q3-3.pdf")
