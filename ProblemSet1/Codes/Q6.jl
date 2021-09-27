using StatsBase
using Plots

Parameters = [[[0 0 ; 0 0.16], [0 ; 0]],
        [[0.85 0.04 ; -0.04 0.85], [0 ; 1.6]],
        [[0.20 -0.26 ; 0.23 0.22], [0 ; 1.6]],
        [[-0.15 0.28 ; 0.26 0.24], [0 ; 0.44]]]

f1(vector, Pars) = Pars[1][1] * vector + Pars[1][2]
f2(vector, Pars) = Pars[2][1] * vector + Pars[2][2]
f3(vector, Pars) = Pars[3][1] * vector + Pars[3][2]
f4(vector, Pars) = Pars[4][1] * vector + Pars[4][2]

function RandomBarnsleyFern(P, PointNums)
    Functions = [f1, f2, f3, f4]
    X = [0.0 for i = 1:PointNums]
    Y = [0.0 for i = 1:PointNums]
    for i in 1:PointNums
        point = [rand() ; rand()]
        for j in 1:P
            s = sample([1,2,3,4], Weights([0.01,0.85,0.07,0.07]))
            point = Functions[s](point, Parameters)
        end
        X[i]= point[1]
        Y[i]= point[2]
    end
    return X, Y
end

x, y = RandomBarnsleyFern(100, 1000000)
scatter(x,y,markersize = 0.0001, legend = false, border=:none, color =:green, dpi=500)
savefig("C:\\Users\\Yaghoub\\Documents\\GitHub\\A-few-fractals-in-Julia\\Fractals-Fig\\RandomBarnsleyFern.png")
