using Plots

Rotation(x, y, θ) = [cos(θ) -sin(θ);sin(θ) cos(θ)] * [x;y]

function KochCurve(; points, order, θ)
    for n in 1:order
        i = 1
        while i < length(points)
            insert!(points,i+1,(points[i+1] - points[i]) /3 + points[i])
            delta = (points[i+2] - points[i]) *2 /3 + points[i] - points[i+1]
            insert!(points, i+2, Rotation(delta[1], delta[2], θ)+points[i+1])
            insert!(points,i+3,(points[i+3] - points[i]) *2 /3 + points[i])
            i += 4
        end
    end
    return points
end

StartPoints = [[[0.0; 0.0],[10.0; 0.0]],
            [[0.0; 0.0],[10.0; 0.0],[5.0; 8.66],[0.0; 0.0]],
            [[0.0; 0.0],[10.0; 0.0],[5.0; 8.66],[0.0; 0.0]],
            [[0.0; 0.0],[10.0; 0.0],[5.0; 8.66],[0.0; 0.0]]
                ]
Adegs = [pi/3, pi/3, -pi/3, 2*pi/3,]

ylimlist = [(0,5), (0,10),(-3,9), (0,10)]
xlimlist = [(0,10), (0,10),(-1,11), (0,10)]

for j in 1:5
    for i in 1:4
        Parameters = Dict(
                    :points => StartPoints[i],
                    :order => j,
                    :θ => Adegs[i])
        points = KochCurve(;Parameters...)
        plot(hcat(points...)[1,:], hcat(points...)[2,:], legend = false, border=:none, fill =(0), ylims = ylimlist[i], xlims = xlimlist[i])
        savefig("C:\\Users\\Yaghoub\\Documents\\GitHub\\ComputationalPhysics-Fall2021\\ProblemSet1\\Figs\\Q1\\KochCurve$(i)O$(j).png")
    end
end
