using ProgressBars, Plots, StatsBase, Statistics, LaTeXStrings, JLD

function RandomWalker2D(TotalSteps, R₀)
    choices = ((1,0),(0,1),(-1,0),(0,-1))
    TotalWalks = rand(choices,TotalSteps)
    TotalR = zeros(Int,TotalSteps,2)
    R = copy(R₀)

    for i ∈ 1:TotalSteps
        R .+= TotalWalks[i]
        TotalR[i,:] = copy(R)
    end
    return TotalR
end

function RGyration(ensemble,RWNumbers)
    TotalR = collect(Tuple(ensemble[i,:]) for i ∈ 1:RWNumbers)
    MassCenter = (mean(TotalR[1]),mean(TotalR[2]))
    return √(sum(pos -> reduce(.+, (pos .- MassCenter).^2), TotalR)/RWNumbers)
end

function ReturnEnsemble(TotalSteps,RWNumbers,R₀)
    Data = zeros(Int,RWNumbers, TotalSteps,2)
    for i ∈ 1:RWNumbers
        Data[i,:,:] = RandomWalker2D(TotalSteps,R₀)
    end

    RgList = zeros(TotalSteps)

    for i ∈ 1:TotalSteps
        RgList[i] = RGyration(Data[:,i,:],RWNumbers)
    end
    return RgList
end

TotalSteps = 40
RWNumbers = 100
R₀ = [0,0]



scatter(ReturnEnsemble(TotalSteps,RWNumbers,R₀))
