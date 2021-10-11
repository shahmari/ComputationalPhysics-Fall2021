using Plots, LaTeXStrings, Statistics

function Deposition(;len, tot_time)
    surf = [0 for i=1:len]
    AvgList = [0.0 for i=1:tot_time]
    for n in 2:tot_time+1
        randsurf = rand(1:len)
        for i in randsurf
            MAX = FindMax(surf, i, len)
            surf[i] = MAX
        end
        AvgList[n-1] = mean(surf)
    end
    return AvgList
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

function FindMax(surface, index_,L_surf)
    i1 , i2 = sides(index_, L_surf)
    maxlen = max(surface[i1],surface[index_] + 1,surface[i2])
    return maxlen
end
function Linear_fit(;Time, AvgList, tot_time)
    A = [hcat(Time[1:tot_time]) reshape(ones(tot_time), tot_time, 1)]
    b = reshape(AvgList[1:tot_time], tot_time, 1)
    line = (A \ b)
    return line
end

iternum = 1000
Parameters = Dict(:len => 300,
                    :tot_time => 50,)
allAvg = [ [0.0 for i in 1:Parameters[:tot_time]] for j = 1:iternum]
meanAvg = [0.0 for i in 1:Parameters[:tot_time]]
vars = [0.0 for i in 1:Parameters[:tot_time]]
for i in 1:iternum
    AvgList = Deposition(;Parameters...)
    allAvg[i] = AvgList
    meanAvg += AvgList
    print("\r$i")
end
meanAvg /= iternum
for i in 1:Parameters[:tot_time]
    vars[i] = std(hcat(allAvg...)[i,:])
end

Time = 0:Parameters[:tot_time]-1

Paraline = Dict(
                :Time => Time[1:end],
                :AvgList => meanAvg[1:end],
                :tot_time => 50)
Line = Linear_fit(;Paraline...)
X = Time[1]:Time[50]
Y = X .* Line[1] .+ Line[2]
Line[2]
plot(X,Y,c=:black,label = L"y = %$(round(Line[1],digits= 3))x + %$(round(Line[2],digits= 3))")
scatter!(Time[1:end], meanAvg[1:end],
    # xlims = (1, Parameters[:tot_time]),
    c = :steelblue,
    xlabel= L"Time",
    ylabel= L"H_{(t)}",
    title= L"~H_{(t)}-Time~\ (L = %$(Parameters[:len]))",
    label = L"Data\ point",
    yerror = vars)
    # legend = nothing)

savefig("../../Figs/Q2/H-t(L=$(Parameters[:len])).pdf")
