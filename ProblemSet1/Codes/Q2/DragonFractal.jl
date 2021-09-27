using Plots

Rotation(R, θ) = [cos(θ) -sin(θ);sin(θ) cos(θ)] * R

function DragonFractal(;order, θ, points)
    θ₁ = θ[1]; θ₂ = θ[2]
    for num in 1:order
        i = 1
        while i < length(points)
            point = (points[i+1] - points[i]) * sin(abs(θ₁))
            point = Rotation(point, θ₁)  + points[i]
            insert!(points,i+1,point)
            i += 2
            point = (points[i+1] - points[i]) * sin(abs(θ₂))
            point = Rotation(point, θ₂)  + points[i]
            insert!(points,i+1,point)
            i += 2
        end
    end
    return points
end


for i in [1,5,10,15]
    Parameters = Dict(:order => i,
                    :θ => [pi/4, -pi/4],
                    :points => [[0.0; 0.0],[1.0; 0.0],[1.0; 1.0]])

    points = DragonFractal(;Parameters...)
    half = Int(floor(length(points)/2))

    plot(hcat(points...)[1,1:half], hcat(points...)[2,1:half], legend = false, color =:blue, border=:none)
    plot!(hcat(points...)[1,half:end], hcat(points...)[2,half:end], legend = false, color =:red, border=:none)
    savefig("C:\\Users\\Yaghoub\\Documents\\GitHub\\ComputationalPhysics-Fall2021\\ProblemSet1\\Figs\\Q2\\Heighway-dragon-O$i.png")

    Parameters = Dict(:order => i,
                    :θ => [-pi/4, -pi/4],
                    :points => [[0.0; 0.0],[1.0; 0.0],[1.0; 1.0]])
    points = DragonFractal(;Parameters...)
    plot(hcat(points...)[1,:], hcat(points...)[2,:], legend = false, border=:none, color =:black, linewidth=0.4)
    savefig("C:\\Users\\Yaghoub\\Documents\\GitHub\\ComputationalPhysics-Fall2021\\ProblemSet1\\Figs\\Q2\\Lévy-C-curve-O$i.png")
end
