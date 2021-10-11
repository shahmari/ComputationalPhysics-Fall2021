using Plots
using Statistics

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

function Deposition_2DG(len, tot_time, rate, color_step)
    surf = [[] for i=1:len + 2]
    for n in 0:tot_time
        randsurf = rand(2:len + 1, rate)
        for i in randsurf
            if i == 2 && length(surf[i]) >= length(surf[i-1]) - 2
                append!(surf[i-1], n%color_step)
            elseif i == len + 1 && length(surf[i]) >= length(surf[i+1]) - 2
                append!(surf[i+1], n%color_step)
            end
            lens_ = Dict(
                length(surf[i-1])=>i-1,
                length(surf[i])=>i,
                length(surf[i+1])=>i+1,
                    )
            minlen = min(length(surf[i-1]),length(surf[i]),length(surf[i+1]))
            append!(surf[lens_[minlen]], n%color_step)
        end
    end
    return Leveling(surf)
end

# theme(:dark)
# plotlyjs()

surf = Deposition_2DG(200, 100, 1000, 10)

heatmap(hcat(surf...), c = :solar, legend = false)
savefig("../../Figs/Q1/Graphic.pdf")
