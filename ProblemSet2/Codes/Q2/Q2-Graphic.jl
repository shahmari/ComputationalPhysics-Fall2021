using Plots, Statistics, LaTeXStrings

function Leveling(surface)
    maxlen = 0
    surf = surface
    for i in surf
        if length(i) > maxlen
            maxlen = length(i)
        end
    end
    for i in surf
        append!(i, [3 for i=1:maxlen - length(i)])
    end
    return surf
end

function Deposition_2DG(;len, tot_time, rate, color_step)
    surf = [[] for i=1:len]
    for n in 0:tot_time
        randsurf = rand(1:len,rate)
        for i in randsurf
            append!(surf[i], n%color_step)
        end
    end
    return hcat(Leveling(surf)...)
end

# theme(:dark)
# gr()

Parameter = Dict( :len => 200,
                    :tot_time => 100,
                        :rate => 10000,
                            :color_step => 10)

surf = Deposition_2DG(;Parameter...)

heatmap(surf, c = :solar, legend = false, title = L"Surface \ length = %$(Parameter[:len]), Total \ number \ of \ particles = 10^{%$(Int(log10(Parameter[:tot_time]*Parameter[:rate])))}")
savefig("C:\\Users\\Yaghoub\\Documents\\GitHub\\ComputationalPhysics-Fall2021\\ProblemSet2\\Figs\\Q2\\RBDGr.png")
