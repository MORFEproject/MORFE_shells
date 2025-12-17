struct SSMParam
    Φ::Array{Int64,1}
    max_order::Int64
    max_orderNA::Int64
    neig::Int64
    Fmodes::Array{Int64,1}
    Fmult::Array{Float64,1}
    Ffreq::Int64
    nzforce::Int64
    omega_mul::Float64
    style::String
    nK::Int64
    nA::Int64
    nrom::Int64
    nz::Int64
    nm::Int64
    nMat::Int64

end

function SSMParamInitShell(MeshInfo::MeshShell)
    Φ = [1,2] # master mode selection
    max_order = 3 # order for autonomous
    max_orderNA = 3 # order for non-autonomous
    neig = 20 # maximum order for eigenvalue solving
    Fmodes = [1] # excitation mode index
    Fmult = (0.5)*[0.1]
    Ffreq = 1
    nzforce = 0 # number of excitation forces (note that the state space expression is used in the program, even if it is 1, it will be decomposed into 2 conjugates)
    omega_mul = 1.0 # a certain multiple frequency parameter
    style = "c" # complex form "c", graphical form "g"
    nK = MeshInfo.NFREEDOF
    nA = 2 * nK
    nrom = nzforce + 2 * length(Φ)
    nz = 2 * length(Φ)
    nm = length(Φ)
    nMat = nA + 2 * length(Φ)
    return SSMParam(Φ, max_order, max_orderNA, neig, Fmodes, Fmult, Ffreq, nzforce, omega_mul, style, nK, nA, nrom, nz, nm, nMat)
end