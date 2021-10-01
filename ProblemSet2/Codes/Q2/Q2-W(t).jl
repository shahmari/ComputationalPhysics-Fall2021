using Plots, LaTeXStrings, Statistics

function Deposition(;len, tot_time, time_steps, rate)
    Time = ceil.(Int, exp.(0:(tot_time)/(time_steps):tot_time))
    surf = [0 for i=1:len]
    VarList = [0.0 for i=1:time_steps]
    for n in 1:time_steps
        randsurf = rand(1:len,(Time[n+1]-Time[n])*rate)
        for i in randsurf
            surf[i] += 1
        end
        VarList[n] = std(surf)
    end
    return VarList
end

function Linear_fit(;len, tot_time, time_steps, rate)
    A = [hcat(log.(Time)) reshape(ones(time_steps), time_steps, 1)]
    b = reshape(log.(meanVar), time_steps, 1)
    line = (A \ b)
    x = 0:tot_time
    y = x .* line[1] .+ line[2]
    return x, y, line
end

iternum = 1000
Parameters = Dict(
                :len => 200,
                :tot_time => 10,
                :time_steps => 20,
                :rate => 10)
allVar = [[0.0 for i in 1:Parameters[:time_steps]] for j = 1:iternum]
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

Time = ceil.(Int, exp.(0:(Parameters[:tot_time])/(Parameters[:time_steps]):Parameters[:tot_time]))[1:end-1]
X, Y, Line = Linear_fit(;Parameters...)

scatter(log.(Time),log.(meanVar),
    xlabel= L"Log\ Time",
    ylabel= L"Log\ W_{(t)}",
    title= L"Log-Log\ Plot\ of\ ~W_{(t)}-Time~, Particles \ number \approx e^{%$(Parameters[:tot_time])}",
    yerror = vars,
    label = L"Data\ point")
plot!(X,Y,c= :black,label = L"y = %$(round(Line[1],digits= 2))x + %$(round(Line[2],digits= 2))")
savefig("C:\\Users\\Yaghoub\\Documents\\GitHub\\ComputationalPhysics-Fall2021\\ProblemSet2\\Figs\\Q2\\RBDW(t).png")
