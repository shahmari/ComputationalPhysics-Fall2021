using Plots

function SierpinskiTriangle(order)
    lasttringles = [[[0.0; 0.0],[10.0; 0.0],[5.0; 8.66],[0.0; 0.0]]]
    for i in 1:order
        appendTringles = []
        for points in lasttringles
            delta = points[1]
            Trgpoints = (points - fill(delta,4)) / 2
            movecor1 = Trgpoints[2]
            movecor2 = Trgpoints[3]
            Trgpoints += fill(delta,4)
            append!(appendTringles, [Trgpoints])
            append!(appendTringles, [Trgpoints + fill(movecor1, 4)])
            append!(appendTringles, [Trgpoints + fill(movecor2, 4)])
        end
        lasttringles = appendTringles
    end
    return lasttringles
end

function PlotingSPT(lasttringles, order)
    backtrg = [[0.0; 0.0],[10.0; 0.0],[5.0; 8.66],[0.0; 0.0]]
    plot(hcat(backtrg...)[1,:], hcat(backtrg...)[2,:], legend = false, border=:none, dpi = 500)
    for points in lasttringles
        plot!(hcat(points...)[1,:], hcat(points...)[2,:], legend = false, fill = (0, 0.5, :green), dpi = 500)
    end
    plot!(dpi = 500)
    savefig("C:\\Users\\Yaghoub\\Documents\\GitHub\\ComputationalPhysics-Fall2021\\ProblemSet1\\Figs\\Q3\\SPT-O$order.png")
end


for i in [1,3,5,7]
    PlotingSPT(SierpinskiTriangle(i),i)
end
