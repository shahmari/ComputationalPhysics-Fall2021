module ExportData

export AnsembleData, SingleData, MineData

cd(dirname(@__FILE__))
include("Q1-MD.jl")
datapath = "../../Data/Q1/"

using ProgressMeter, JLD, Statistics

function AnsembleData(sn::Integer, rn::Integer, Parameters::Dict)
    stepnum = sn
    runnum = rn
    Prog = Progress(stepnum * runnum)

    TotUColl = zeros(stepnum, runnum)
    TotKColl = zeros(stepnum, runnum)
    TotTColl = zeros(stepnum, runnum)
    TotPColl = zeros(stepnum, runnum)
    TotrsColl = zeros(stepnum, runnum)
    TotlsColl = zeros(stepnum, runnum)

    for j ∈ 1:runnum
        sys = MDSim.init(; Parameters...)
        for i ∈ 1:stepnum
            MDSim.simulate!(sys)
            lsN = count(<=(sys.l / 2), sys.r[1, :]) / sys.N
            rsN = 1 - lsN
            TotlsColl[i, j] = lsN
            TotrsColl[i, j] = rsN
            TotUColl[i, j] = sys.U
            TotKColl[i, j] = sys.K
            TotTColl[i, j] = sys.T
            TotPColl[i, j] = sys.P
            next!(Prog)
            update!(Prog)
        end
    end
    save(datapath * "DataTD.jld", "TotUColl", TotUColl, "TotKColl", TotKColl,
        "TotTColl", TotTColl, "TotPColl", TotPColl, "TotlsColl", TotlsColl, "TotrsColl", TotrsColl)
end

function SingleData(sn::Integer, Parameters::Dict)
    sys = MDSim.init(; Parameters...)

    stepnum = sn
    Prog = Progress(stepnum)
    PosColl = zeros(stepnum, 2, Parameters[:N])
    VelColl = zeros(stepnum, 2, Parameters[:N])

    UColl = zeros(stepnum)
    KColl = zeros(stepnum)
    TColl = zeros(stepnum)
    PColl = zeros(stepnum)

    for i ∈ 1:stepnum
        MDSim.simulate!(sys)
        PosColl[i, :, :] = sys.r
        VelColl[i, :, :] = sys.v
        UColl[i] = sys.U
        KColl[i] = sys.K
        TColl[i] = sys.T
        PColl[i] = sys.P
        next!(Prog)
        update!(Prog)
    end
    save(datapath * "DataSD.jld",
        "PosColl", PosColl, "VelColl", VelColl,
        "UColl", UColl, "KColl", KColl, "TColl", TColl, "PColl", PColl)
end

function MineData(TotUColl::Matrix{T}, TotKColl::Matrix{T}, TotTColl::Matrix{T}, TotPColl::Matrix{T}, TotlsColl::Matrix{T}, TotrsColl::Matrix{T}) where {T<:AbstractFloat}
    TUC = [[TotUColl[i, :]...] for i ∈ 1:size(TotUColl)[1]]
    for i ∈ TUC
        deleteat!(i, findall(x -> isnan(x), i))
    end
    TKC = [[TotKColl[i, :]...] for i ∈ 1:size(TotKColl)[1]]
    for i ∈ TKC
        deleteat!(i, findall(x -> isnan(x), i))
    end
    TTC = [[TotTColl[i, :]...] for i ∈ 1:size(TotTColl)[1]]
    for i ∈ TTC
        deleteat!(i, findall(x -> isnan(x), i))
    end
    TPC = [[TotPColl[i, :]...] for i ∈ 1:size(TotPColl)[1]]
    for i ∈ TPC
        deleteat!(i, findall(x -> isnan(x), i))
    end
    TRSC = [[TotrsColl[i, :]...] for i ∈ 1:size(TotrsColl)[1]]
    for i ∈ TRSC
        deleteat!(i, findall(x -> isnan(x), i))
    end
    TLSC = [[TotlsColl[i, :]...] for i ∈ 1:size(TotlsColl)[1]]
    for i ∈ TLSC
        deleteat!(i, findall(x -> isnan(x), i))
    end

    save(datapath * "MeanData.jld", "MEANUColl", mean.(TUC), "MEANKColl", mean.(TKC),
        "MEANTColl", mean.(TTC), "MEANPColl", mean.(TPC), "MEANlsColl", mean.(TLSC), "MEANrsColl", mean.(TRSC))
    save(datapath * "STDData.jld", "STDUColl", std.(TUC), "STDKColl", std.(TKC),
        "STDTColl", std.(TTC), "STDPColl", std.(TPC), "STDlsColl", std.(TLSC), "STDrsColl", std.(TRSC))
end

end