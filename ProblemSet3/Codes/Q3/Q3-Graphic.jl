using Plots, LaTeXStrings, Statistics

function Deposition_2DG(; len, tot_time, dep_rate, color_step)
    surface = InitialSurf(len)
    for n in 0:tot_time
        randsurf = rand(1:len, dep_rate)
        for i in randsurf
            surface = check_and_depose(surface, i, len, (n%color_step) + 1)
        end
    end
    return surface
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

Parameters = Dict(:len => 300,
                    :tot_time => 1000,
                        :dep_rate => 50,
                            :color_step => 100)

surf = Deposition_2DG(;Parameters...)

heatmap(transpose(hcat(surf...)), c = :solar, legend = false)
savefig("../../Figs/Q3/Graphic.pdf")
