module ExportData

export AnsembleData, SingleData

cd(dirname(@__FILE__))
include("Q1-MD.jl")

using Statistics: mean
using ProgressMeter

function AnsembleData(sn::Integer, rn::Integer, Parameters::Dict)
    stepnum = sn
    runnum = rn
    Prog = Progress(stepnum * runnum)

    TotUColl = zeros(stepnum, runnum)
    TotKColl = zeros(stepnum, runnum)
    TotTColl = zeros(stepnum, runnum)
    TotPColl = zeros(stepnum, runnum)

    for j ∈ 1:runnum
        sys = MDSim.init(; Parameters...)
        for i ∈ 1:stepnum
            MDSim.simulate!(sys)
            TotUColl[i, j] = sys.U
            TotKColl[i, j] = sys.K
            TotTColl[i, j] = sys.T
            TotPColl[i, j] = sys.P
            next!(Prog)
            update!(Prog)
        end
    end
    return mean(TotUColl, dims = 2), mean(TotKColl, dims = 2), mean(TotTColl, dims = 2), mean(TotPColl, dims = 2)
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
    return UColl, KColl, TColl, PColl, PosColl, VelColl
end

end