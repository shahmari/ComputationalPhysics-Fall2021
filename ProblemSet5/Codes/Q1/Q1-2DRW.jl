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

histogram2d(collect(Tuple(Data[end,:,i]) for i ∈ 1:1000), nbins = 20,framestyle = :box, background=:steelblue)
scatter!(collect(Tuple(Data[i,:]) for i ∈ 1:TotalSteps))

Data = zeros(1000,1000,2)
Data[i,:,:] for i = 1:1000 .= RandomWalker2D(1000, R₀)
Data = cat([RandomWalker2D(10000, R₀) for i = 1:1000]...; dims = 3)
Data[end,:,:]
