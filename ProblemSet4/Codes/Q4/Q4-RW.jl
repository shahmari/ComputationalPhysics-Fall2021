using ProgressBars, Plots, StatsBase, Statistics, LaTeXStrings, JLD


function RandomWalker(P,TimeSteps)
    X = 0
    for i in 1:TimeSteps
        X += sample([1,-1], Weights([P,1-P]))
    end
    return X
end

function Linear_fit(X, Y)
    A = [hcat(X) reshape(ones(length(X)), length(X), 1)]
    b = reshape(hcat(Y), length(Y), 1)
    line = (A \ b)
    return line
end

PList = [0.1,0.3,0.5,0.7,0.9]
N = 100000
TimeSteps = 10000
RES = [[] for i = 1:5]

for n in 1:5
    P = PList[n]
    for i in ProgressBar(1:N)
        push!(RES[n],RandomWalker(P,TimeSteps))
    end
    push!(Avg,mean(RES[n]))
end

Line = round.(Linear_fit(PList,Avg))

p1 = histogram(RES[1],bins = 50, label = L"P = %$(PList[1]), \ l = \tau = 1", xlabel=L"\langle X_{final} \rangle_{(P)}", ylabel = L"Number", c = :steelblue)
p2 = histogram(RES[2],bins = 50, label = L"P = %$(PList[2]), \ l = \tau = 1", xlabel=L"\langle X_{final} \rangle_{(P)}", ylabel = L"Number", c = :purple)
p3 = histogram(RES[3],bins = 50, label = L"P = %$(PList[3]), \ l = \tau = 1", xlabel=L"\langle X_{final} \rangle_{(P)}", ylabel = L"Number", c = :red)
p4 = histogram(RES[4],bins = 50, label = L"P = %$(PList[4]), \ l = \tau = 1", xlabel=L"\langle X_{final} \rangle_{(P)}", ylabel = L"Number", c = :gold)
p5 = histogram(RES[5],bins = 50, label = L"P = %$(PList[5]), \ l = \tau = 1", xlabel=L"\langle X_{final} \rangle_{(P)}", ylabel = L"Number", c = :green)
p6 = begin
    plot(PList,Avg,label = L"Y = %$(round(Line[1],digits=3)) X %$(round(Line[2],digits=3))", c = :black)
    scatter!(PList,Avg,label = L"Data\ Points",legend = 140,xlabel = L"P",ylabel = L"\langle \overline{X_{final}} \rangle_{(P)}",c = :steelblue)
end

Pmain1 = begin
    plot(plot(plot(p1,p2,layout=(1,2)),plot(p3,p4,layout=(1,2)),plot(p5,p6,layout=(1,2)),layout=(3,1)))
    plot!(plot_title = L"Distribution\ of\ X_{final}\ (N = %$N\ , Time\ Steps=%$TimeSteps)",size = (1000,1000))
end
savefig(Pmain1,"../../Figs/Q4/Q4-tot.pdf")
