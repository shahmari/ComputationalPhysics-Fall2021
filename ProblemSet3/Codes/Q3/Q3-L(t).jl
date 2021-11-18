using Plots, LaTeXStrings, Statistics, ProgressBars, JLD
cd(dirname(@__FILE__))

function Deposition(; len, tot_time, dep_rate, color_step)
    surface = InitialSurf(len)
    AllLens = []
    for n in 0:tot_time
        randsurf = rand(1:len, dep_rate)
        for i in randsurf
            surface = check_and_depose(surface, i, len, (n%color_step) + 1)
        end
        append!(AllLens, Check_Dlen(surface))
    end
    return AllLens
end

function Check_Dlen(surf)
    Dlen = []
    for row in surf
        yindex = []
        for num in 1:length(row)
            if row[num] != 0
                append!(yindex, num)
            end
        end
        append!(Dlen, max(yindex...) - min(yindex...))
    end
    return max(Dlen...)
end

function check_and_depose(surf, index, len, addval) #checking down and sides and first layer
    i1 , i2 = sides(index, len)
    if surf[end][index] != 0
        newsurf = [0 for j = 1:len]
        newsurf[index] = addval
        push!(surf, newsurf)
    end
    for num in 0:length(surf) -1
        if surf[end - num][index] == 0 && (surf[end - num][i1] != 0 || surf[end - num][i2] != 0)
            surf[end - num][index] = addval
            break
        elseif length(surf) > num + 1 && surf[end - 1 - num][index] != 0
            surf[end - num][index] = addval
            break
        end
    end
    return surf
end

function sides(n, L)
    if n == L
        return n-1 , 1
    elseif n == 1
        return L , n+1
    else
        return n-1, n+1
    end
end

function InitialSurf(len)
    pos1 = [0 for i=1:len]
    pos1[floor(Int, len/2)] = 1
    return [pos1]
end

function Linear_fit(;Time, AvgList, steps)
    A = [hcat(Time) reshape(ones(steps), steps, 1)]
    b = reshape(AvgList, steps, 1)
    line = (A \ b)
    return line
end

Parameters = Dict(:len => 200,
                    :tot_time => 403,
                        :dep_rate => 100,
                            :color_step => 100)

up = 5.35; low = 1; steps = 100; SelctedSpan = 56:95; SelectedSteps = length(SelctedSpan);
Time = round.(Int,exp.(low:(up-low)/(steps-1):up))[SelctedSpan]
RawTime = hcat(low:(up-low)/(steps-1):up)[SelctedSpan]

iternum = 1000
allVal = [[0.0 for i ∈ 1:SelectedSteps] for j ∈ 1:iternum]
meanVal = [0.0 for i ∈ 1:SelectedSteps]
vars = [0.0 for i ∈ 1:SelectedSteps]
for i ∈ ProgressBar(1:iternum)
    VarList = Deposition(;Parameters...)[Time]
    allVal[i] = VarList
    meanVal += VarList
end
meanVal /= iternum
for i ∈ 1:SelectedSteps
    vars[i] = std(hcat(allVal...)[i,:])
end

Paraline = Dict(:Time => log.(Time),
                :AvgList => log.(meanVal),
                :steps => SelectedSteps)
Line = Linear_fit(;Paraline...)
X = RawTime
Y = X .* Line[1] .+ Line[2]

plot(X,Y,c = :black,label = L"y = %$(round(Line[1],digits= 3))x + %$(round(Line[2],digits= 3))")
scatter!(RawTime, log.(meanVal),
    c = :steelblue,
    xlabel= L"\mathrm{ln}(Time)",
    ylabel= L"\mathrm{ln}({L_C}_{(t)})",
    title= L"Log-Log\ Plot\ of\ L_C\ (L = %$(Parameters[:len]),\ %$iternum\ runs)",
    label = L"Data\ point",
    ribbon = log.(vars),
    yerror = log.(vars),
    legend = 150)

savefig("../../Figs/Q3/Q3-L(t).pdf")
save("../../Data/Q3/Q3-TotData.jld","Data", allVal)