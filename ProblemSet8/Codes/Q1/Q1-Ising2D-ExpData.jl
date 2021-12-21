using ProgressMeter, JLD

include("Q1-Ising2D.jl")

cd(dirname(@__FILE__))

datapath = "../../Data/Q1/"

rnum = 30
Llist = hcat(4:20)
βzoom = hcat(0.35:0.002:0.45)
EXPData = zeros(4, length(βzoom), length(Llist), rnum)
PRG = Progress(length(βzoom) * rnum)
for n ∈ 1:rnum
    for i ∈ 1:length(βzoom)
        for l ∈ 1:length(Llist)
            L = Llist[l]
            β = βzoom[i]
            Parameters = Dict(
                :dim => L,
                :MCLSize => 1000,
                :NLSize => 1000,
                :SkipNum => 200,
                :ProgBar => false)
            EXPData[:, i, l, n] = [Ising2D.IsingModel(β; Parameters...)...]
        end
        next!(PRG)
    end
end
save(datapath * "EXPData.jld", "EXPData", EXPData, "βzoom", βzoom, "Llist", Llist)