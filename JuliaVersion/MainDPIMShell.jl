using Pkg

Pkg.activate(".")
Pkg.instantiate()
using DelimitedFiles
using SparseArrays
using ExtendableSparse
using LinearAlgebra
using Arpack
using BenchmarkTools
using Combinatorics
using MAT
using Printf

include("./FEM/ReadMesh.jl")
include("./DPIM/DPIMParam.jl")
include("./FEM/MaterialParam.jl")
include("./DPIM/ParamInit.jl")
include("./FEM/ElementIntepolation.jl")
include("./FEM/Assemble.jl")
include("./FEM/NonlinearForce.jl")
include("./DPIM/dpim.jl")
include("./DPIM/dpim_routines.jl")
include("./DPIM/realification.jl")
include("./Tool/export_solution.jl")

# FEM mesh document
NodeFileName = "./NumericalExample/ShallowShell_M/NLIST.lis"
ElementFileName = "./NumericalExample/ShallowShell_M/ELIST.lis"
# add thickness column to the node file
addThicknessColumn(NodeFileName, "TestThickness.txt", 0.01, 0.01)
NodeFileName = "TestThickness.txt"
# read fem information
MeshInfo = ReadShell(NodeFileName, ElementFileName)
# set DPIM parameters
DPIMParams = DPIMParamInitShell(MeshInfo)
# set force parameters
ForceVectorFull = zeros(Float64, MeshInfo.NDOF)
MeshInfo.ForceNodeDOF = [1]
ForceVectorFull[MeshInfo.ForceNodeDOF] = DPIMParams.Fmult
MeshInfo.ForceVector = ForceVectorFull[MeshInfo.FREEDOF]
# set material parameters
MatParams = MatParamInitShell(MeshInfo.NodeCoord)
# precalculate matrix based on total lagrange
PreCalculate!(MeshInfo, MatParams)
# parameterize spectral subspace manifold
result = @time begin
    ϕ, MF, Cp = dpim(MeshInfo, DPIMParams, MatParams)
    realification!(Cp,DPIMParams)
    howmany=count_terms_dyn(Cp,DPIMParams)
    # store dynamic and mapping information
    mappings,mappings_vel,mappings_modal,mappings_modal_vel,Avector,fdyn = store_dyn_and_map(Cp,DPIMParams,howmany, MeshInfo, MatParams)
    # write reduced dynamic information
    write_rdyn(DPIMParams,Cp) 
end
# save data that can be used in matlab code
SaveData(mappings,mappings_vel,mappings_modal,mappings_modal_vel,Avector,fdyn,howmany,MeshInfo.ForceNode)
println("done!")
