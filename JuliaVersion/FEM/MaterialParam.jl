struct MaterialParamShell
    E::Float64
    ν::Float64
    ρ::Float64
    α::Float64
    β::Float64
    t::Array{Float64,1}
end

function MatParamInitShell(NodeCoord)
    E = 70e9
    ν = 0.33
    ρ = 2700
    # α = 1*0.001*133.7945
    α = 0
    β = 0
    t = NodeCoord[:,end]
    return MaterialParamShell(E, ν, ρ, α, β, t)
end