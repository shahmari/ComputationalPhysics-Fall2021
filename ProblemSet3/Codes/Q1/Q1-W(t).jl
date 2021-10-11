using Plots, LaTeXStrings, Statistics

function Deposition(;len, tot_time, time_steps)
    Time = ceil.(Int, exp.(2:(tot_time-2)/(time_steps):tot_time))
    surf = [0 for i=1:len]
    VarList = [0.0 for i=1:time_steps]
    for n in 2:time_steps+1
        randsurf = rand(1:len,(Time[n]-Time[n-1]))
        for i in randsurf
            index = FindLeast(surf, i, len)
            surf[index] += 1
        end
        VarList[n-1] = std(surf)
    end
    return VarList
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

function FindLeast(surface, index_,L_surf)
    i1 , i2 = sides(index_, L_surf)
    lens = Dict(surface[i1]=>i1,
        surface[index_]=>index_,
        surface[i2]=>i2)
    minlen = min(surface[i1],surface[index_],surface[i2])
    return lens[minlen]
end

function Linear_fit(;Time, VarList, time_steps)
    A = [hcat(log.(Time[1:time_steps])) reshape(ones(time_steps), time_steps, 1)]
    b = reshape(log.(VarList[1:time_steps]), time_steps, 1)
    line = (A \ b)
    return line
end

# function FindLine(tot_steps, Time, VarList)
#     retline = [0.0 0.0]
#     for i in 5:tot_steps
#         Paraline = Dict(
#                         :Time => Time,
#                         :VarList => VarList,
#                         :time_steps => i
#                             )
#         Line = Linear_fit(;Paraline...)
#         absl = abs(Line[1]*Time[i] + Line[2] - VarList[i])
#         standev = std(Line[1].*Time .+ Line[2] .- VarList)
#         if absl > 0.2 && standev > 1 && retline == [0.0 0.0]
#             retline = Line
#         end
#         if absl > 1 && standev > 3
#             return retline, i
#         end
#     end
# end

# theme(:dark)
# gr()


iternum = 1000
Parameters = Dict(:len => 300,
                    :tot_time => 15,
                        :time_steps => 50)
allVar = [ [0.0 for i in 1:Parameters[:time_steps]] for j = 1:iternum]
meanVar = [0.0 for i in 1:Parameters[:time_steps]]
vars = [0.0 for i in 1:Parameters[:time_steps]]
for i in 1:iternum
    VarList = Deposition(;Parameters...)
    allVar[i] = VarList
    meanVar += VarList
    print("\r$i")
end
meanVar /= iternum
for i in 1:Parameters[:time_steps]
    vars[i] = std(log.(hcat(allVar...))[i,:])
end

Time = exp.(0:(Parameters[:tot_time])/(Parameters[:time_steps]-1):Parameters[:tot_time])
scatter(log.(Time[10:end]), log.(meanVar[10:end]),
    # xlims = (1, Parameters[:tot_time]),
    xlabel= L"Log\ Time",
    ylabel= L"Log\ W_{(t)}",
    title= L"Log-Log\ Plot\ ~W_{(t)}-Time~\ (L = %$(Parameters[:len]))",
    label = L"Data\ point",
    yerror = vars)
    # legend = nothing)

Paraline = Dict(
                :Time => Time[10:end],
                :VarList => meanVar[10:end],
                :time_steps => 30)
Line = Linear_fit(;Paraline...)
X = log.(Time[10]:Time[50])
Y = X .* Line[1] .+ Line[2]

plot!(X,Y,c=:black,label = L"y = %$(round(Line[1],digits= 2))x + %$(round(Line[2],digits= 2))")

savefig("../../Figs/Q1/W-t(L=$(Parameters[:len])).pdf")
