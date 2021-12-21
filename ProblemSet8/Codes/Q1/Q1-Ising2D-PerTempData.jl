using ProgressMeter, JLD

include("Q1-Ising2D.jl")

cd(dirname(@__FILE__))

datapath = "../../Data/Q1/"

βList = [hcat(0.0:0.02:0.25); hcat(0.26:0.0075:0.55); hcat(0.56:0.02:1.0)]
Llist = [5, 10, 15]
runnum = 100
Data = zeros(4, length(βList), runnum, length(Llist))

PRG = Progress(length(βList) * runnum)
for n ∈ 1:runnum
    for i ∈ 1:length(βList)
        for l ∈ 1:length(Llist)
            L = Llist[l]
            β = βList[i]
            Parameters = Dict(
                :dim => L,
                :MCLSize => 1500,
                :NLSize => 1000,
                :SkipNum => 500,
                :ProgBar => false)
            Data[:, i, n, l] = [Ising2D.IsingModel(β; Parameters...)...]
        end
        next!(PRG)
    end
end
save(datapath * "PTData.jld", "Data", Data, "βList", βList)
