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

order = 2^11
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
scatter(x,y,markersize = 0.000001, legend = false, border=:none, dpi=300)
savefig("C:\\Users\\Yaghoub\\Documents\\GitHub\\A-few-fractals-in-Julia\\Fractals-Fig\\SKT.png")
