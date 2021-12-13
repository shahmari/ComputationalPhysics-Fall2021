module Ising2D

export NetworkDynamic, IsingModel

using Statistics, ProgressMeter

# function FixBounds(dim::Integer, indices::Vector{Vector{Int}})
#     for i ∈ 1:4
#         for j ∈ 1:2
#             if indices[i][j] == dim + 1
#                 indices[i][j] = 1
#             elseif indices[i][j] == 0
#                 indices[i][j] = dim
#             end
#         end
#     end
#     return indices
# end

function FindEnergy(Network::Matrix{Int8}, dim::Integer)
    E = 0.0
    for i ∈ 1:dim
        for j ∈ 1:dim
            Neighbors = [[1 + i % dim, j], [i, 1 + j % dim], [(-2 + i + dim) % dim + 1, j], [i, (-2 + j + dim) % dim + 1]]
            E += -sum(Network[i, j] * Network[CartesianIndex.(Tuple.(Neighbors))]) / 2
        end
    end
    return E / dim ^2
end

function FindMagnetization(Network::Matrix{Int8})
    return mean(Network)
end

function InitialCondition(dim::Integer)
    return rand([Int8(-1), Int8(1)], dim, dim)
end

function NetworkDynamic(β::Real; dim::Integer = 150, MCLSize::Integer = 200, NLSize::Integer = 5000, ProgBar::Bool = true)
    if ProgBar == true
        Prog = Progress(MCLSize * NLSize)
    end
    Networks = Array{Matrix{Int8},1}(undef, MCLSize)
    Network = InitialCondition(dim)
    for MCL ∈ 1:MCLSize
        for NL ∈ 1:NLSize
            if ProgBar == true
                next!(Prog)
            end
            i, j = rand(1:dim, 2)
            Neighbors = [
                [1 + i % dim, j],
                [i, 1 + j % dim],
                [(-2 + i + dim) % dim + 1, j],
                [i, (-2 + j + dim) % dim + 1]
            ]
            ΔE = 2 * Network[i, j] * sum(Network[CartesianIndex.(Tuple.(Neighbors))])
            if ΔE <= 0
                Network[i, j] *= -1
            elseif rand() < exp(-ΔE * β)
                Network[i, j] *= -1
            end
        end
        Networks[MCL] = copy(Network)
    end
    return Networks
end

function MonteCarloStep(StepSize::Integer, Network::Matrix{Int8}, dim::Integer, β::Real)
    ΔEₙₑₜ = 0
    ΔMₙₑₜ = 0
    for NL ∈ 1:StepSize
        i, j = rand(1:dim, 2)
        Neighbors = [
            [1 + i % dim, j],
            [i, 1 + j % dim],
            [(-2 + i + dim) % dim + 1, j],
            [i, (-2 + j + dim) % dim + 1]
        ]
        ΔE = 2 * Network[i, j] * sum(Network[CartesianIndex.(Tuple.(Neighbors))])
        if ΔE <= 0
            Network[i, j] *= -1
            ΔEₙₑₜ += ΔE / dim^2
            ΔMₙₑₜ += 2 * Network[i, j] / dim^2
        elseif rand() < exp(-ΔE * β)
            Network[i, j] *= -1
            ΔEₙₑₜ += ΔE / dim^2
            ΔMₙₑₜ += 2 * Network[i, j] / dim^2
        end
    end
    return ΔEₙₑₜ, ΔMₙₑₜ
end


function IsingModel(β::Real; dim::Integer = 20, MCLSize::Integer = 1000000, NLSize::Integer = 10, SkipNum::Integer = 500000, ProgBar::Bool = false)
    if ProgBar == true
        Prog = Progress(MCLSize)
    end
    Network = InitialCondition(dim)
    E = FindEnergy(Network, dim)
    M = FindMagnetization(Network)
    EₙₑₜList = Array{Float64,1}(undef, MCLSize)
    MₙₑₜList = Array{Float64,1}(undef, MCLSize)
    for MCL ∈ 1:MCLSize
        if ProgBar == true
            next!(Prog)
        end
        ΔEₙₑₜ, ΔMₙₑₜ = MonteCarloStep(NLSize, Network, dim, β)
        E += ΔEₙₑₜ
        M += ΔMₙₑₜ
        EₙₑₜList[MCL] = E
        MₙₑₜList[MCL] = abs(M)
    end
    Ē, M̄ = mean(EₙₑₜList[SkipNum:end]), mean(MₙₑₜList[SkipNum:end])
    Cᵥ, 𝑋 = β^2 * var(EₙₑₜList[SkipNum:end]), β * var(MₙₑₜList[SkipNum:end])
    return Ē, M̄, Cᵥ, 𝑋
end

end