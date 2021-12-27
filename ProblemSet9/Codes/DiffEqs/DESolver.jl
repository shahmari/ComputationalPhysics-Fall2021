module DiffEqSolver

export EulerDES

function EulerDES(; ḟ::Function, y₀::Real, x₀::Float64, x₁::Float64, h::Float64)
    xcoll = hcat(x₀:h:x₁)
    ycoll = Float64[y₀]
    y = y₀
    for x ∈ xcoll
        y += h * ḟ(x, y)
        push!(ycoll, y)
    end
    return xcoll, ycoll
end

end


ts = collect(range(0.0, 5.0, step = 0.01))



Vector{Float64}(undef, length(hcat(0.0:0.01:5.0)))