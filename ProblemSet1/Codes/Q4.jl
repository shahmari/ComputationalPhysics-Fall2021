using Plots

func1(R) = R/2
func2(R) = R/2 + [1 ; 0]
func3(R) = R/2 + [0 ; 1]
P = 200
num = 500000

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

scatter(x,y,markersize = 0.000001, legend = false, border=:none, dpi=500)
savefig("C:\\Users\\Yaghoub\\Documents\\GitHub\\A-few-fractals-in-Julia\\Fractals-Fig\\SRT.png")
