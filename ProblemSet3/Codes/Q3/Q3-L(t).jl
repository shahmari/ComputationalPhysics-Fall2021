using Plots, LaTeXStrings, Statistics

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

function Linear_fit(;Time, AvgList, tot_time)
    A = [hcat(Time[1:tot_time]) reshape(ones(tot_time), tot_time, 1)]
    b = reshape(AvgList[1:tot_time], tot_time, 1)
    line = (A \ b)
    return line
end

Parameters = Dict(:len => 300,
                    :tot_time => 50,
                        :dep_rate => 1000,
                            :color_step => 100)


iternum = 200
allVal = [ [0.0 for i in 0:Parameters[:tot_time]] for j = 1:iternum]
meanVal = [0.0 for i in 0:Parameters[:tot_time]]
vars = [0.0 for i in 0:Parameters[:tot_time]]
for i in 1:iternum
    VarList = Deposition(;Parameters...)
    allVal[i] = VarList
    meanVal += VarList
    print("\r$i")
end
meanVal /= iternum
for i in 1:Parameters[:tot_time] + 1
    vars[i] = std(hcat(allVal...)[i,:])
end

Paraline = Dict(
                :Time => hcat(0:Parameters[:tot_time]),
                :AvgList => meanVal,
                :tot_time => 30)
Line = Linear_fit(;Paraline...)
X = 0:Parameters[:tot_time]
Y = X .* Line[1] .+ Line[2]

plot(X,Y,c = :black,label = L"y = %$(round(Line[1],digits= 3))x + %$(round(Line[2],digits= 3))")
scatter!(X, meanVal,
    c = :steelblue,
    xlabel= L"Time",
    ylabel= L"Transverse\ width",
    title= L"Transverse\ width-Time~\ (L = %$(Parameters[:len]))",
    label = L"Data\ point",
    yerror = vars,
    legend = 150)

savefig("../../Figs/Q3/Q3-L(t).pdf")
