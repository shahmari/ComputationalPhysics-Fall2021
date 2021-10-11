using Plots, LaTeXStrings, Statistics

function Deposition_2DG(; len, tot_time, dep_rate, color_step)
    surf = [[] for i=1:len]
    for i in surf
        push!(i, [0, 0])
    end
    for n in 0:tot_time
        randsurf = rand(1:len, dep_rate)
        for i in randsurf
            maxlen = max(surf[sides(i, len)[1]][end][1], surf[sides(i, len)[2]][end][1], surf[i][end][1]+1)
            push!(surf[i], [maxlen, n%color_step])
        end
    end
    return Leveling(surf, len)
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

function Leveling(surf, len)
    lasth = []
    for i in surf
        append!(lasth, i[end][1])
    end

    surface = zeros(Int,(len, max(lasth...)+1))

    for i in 1:len
        for j in surf[i]
            surface[i, j[1]+1] = j[2]
        end
    end
    return transpose(hcat(surface))
end

Parameters = Dict(:len => 200,
                    :tot_time => 100,
                        :dep_rate => 1000,
                            :color_step => 10)

surf = Deposition_2DG(; Parameters...)

# theme(:dark)
# gr()

heatmap(surf, c = :solar, legend = false)
savefig("../../Figs/Q2/Graphic.pdf")
