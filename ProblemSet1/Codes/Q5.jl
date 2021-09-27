using Plots

function khayyam_triangle(lentingle)
    output = []
    for i in 1:lentingle
        push!(output, [])
        for j in 1:i
            if j == 1
                append!(output[i], 1)
            elseif j == i
                append!(output[i], 1)
            else j != 1 && j != i
                add = output[i - 1][j] + output[i - 1][j - 1]
                append!(output[i], add)
            end
        end
    end
    return output
end

function SPT_KhPT(order)
    order = 2^order
    triangle = khayyam_triangle(order)
    x = []
    y = []
    for i in 1:order
        for j in 1:i
            if triangle[i][j]%2 != 0
                push!(x, j - i/2)
                push!(y, -i)
            end
        end
    end
    return x, y
end

msize = [7,1,0.1]
olist = [5, 7, 10]
for i in 1:3
    x,y = SPT_KhPT(olist[i])
    scatter(x,y,markersize = msize[i], legend = false, border=:none, dpi=300)
    savefig("C:\\Users\\Yaghoub\\Documents\\GitHub\\ComputationalPhysics-Fall2021\\ProblemSet1\\Figs\\Q5\\SPTKP-O$(i).png")
end
