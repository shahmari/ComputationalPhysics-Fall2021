using ProgressBars, Plots, StatsBase, Statistics, LaTeXStrings, JLD

choices = ((1,0),(0,1),(-1,0),(0,-1))

#Diffusion-limited aggregation
function Deposing(;Network,RWpos,cap,ParticlesNumber,H1_0,H1,H2)
    Network[RWpos...,2] += 1
    if Network[RWpos...,2] == cap
        Network[RWpos...,1] = ParticlesNumber
        H1 = max(RWpos[2] + H1_0, H1)
        H2 = H1 + 5
    end
    check = false
    return Network, H1, H2, check
end

function NewRWpos(RWpos,L)
    retRWpos = RWpos .+ rand(choices)
    if retRWpos[1] == L +1
        retRWpos[1] = 1
    elseif retRWpos[1] == 0
        retRWpos[1] = L
    end
    return retRWpos
end

function DLASimulation(L, H1_0, cap)
    Network = zeros(L,floor(Int,L/2),2)
    Network[:,1,1] = ones(L)
    H2_0 = H1_0 + 5
    H1 = H1_0
    H2 = H2_0
    RWpos = [rand(1:L), H1]
    ParticlesNumber = 0
    while H2 < floor(Int,L/2)
        ParticlesNumber += 1
        print("\r$ParticlesNumber   $H2     ")
        RWpos = [rand(1:L), H1]
        check = true
        while RWpos[2] < H2 && check
            RWpos = NewRWpos(RWpos,L)
            if RWpos[2] <= H1 - H1_0 + 2
                Parameters = Dict(:Network => Network, :RWpos => RWpos,
                                    :cap => cap, :ParticlesNumber => ParticlesNumber,
                                        :H1_0 => H1_0, :H1 => H1, :H2 => H2)
                for neighbor in choices
                    if RWpos[1] + neighbor[1] == 0
                        if Network[L,RWpos[2] + neighbor[2],1] != 0
                            Network, H1, H2, check = Deposing(;Parameters...)
                            break
                        end
                    elseif RWpos[1] + neighbor[1] == L + 1
                        if Network[1,RWpos[2] + neighbor[2],1] != 0
                            Network, H1, H2, check = Deposing(;Parameters...)
                            break
                        end
                    elseif Network[RWpos .+ neighbor ...,1] != 0
                        Network, H1, H2, check = Deposing(;Parameters...)
                        break
                    end
                end
            end
        end
    end
    return Network, ParticlesNumber
end

L = 200
H1_0 = 10
cap = 3

Network, ParticlesNumber = DLASimulation(L, H1_0, cap)
save("../../Data/Q2/Q2-HM.jld", "Network", Network)
heatmap(rotr90(Network[:,:,1]), c = cgrad(:turbo, rev = false), legend = nothing,
    title = L"Distribution\ of\ %$ParticlesNumber\ RW\ in\ DLA")
savefig("../../Figs/Q2/Q2-HM.pdf")
