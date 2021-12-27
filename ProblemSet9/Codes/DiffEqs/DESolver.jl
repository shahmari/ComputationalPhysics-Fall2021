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