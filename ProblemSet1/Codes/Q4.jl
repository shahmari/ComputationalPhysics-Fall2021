using Plots

func1(R) = R/2
func2(R) = R/2 + [1 ; 0]
func3(R) = R/2 + [0 ; 1]

function RandomSPT(P, num)
    x = []
    y = []
    for i in 1:num
        point = [rand() ; rand()]
        for j in 1:P
            point = rand([func1(point), func2(point), func3(point)])
        end
        push!(x, point[1])
        push!(y, point[2])
    end
    return x, y
end

for i in [5, 500]
    for j in [10000,100000]
        x,y = RandomSPT(i, j)

        scatter(x,y,markersize = 10^4/j, legend = false, border=:none, dpi=200, label = "Points number = $j Operations number = $i")
        savefig("C:\\Users\\Yaghoub\\Documents\\GitHub\\ComputationalPhysics-Fall2021\\ProblemSet1\\Figs\\Q4\\RSPTp$(i)num$(j).png")
    end
end
