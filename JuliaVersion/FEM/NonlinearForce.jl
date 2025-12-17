function fillrhsG!(MeshInfo::MeshShell,MatParams::MaterialParamShell,Cp::Vector{Parametrisation},p::Int64)
    for p1 in 1:p-1, p2 in 1:p-1
        if (p1+p2)==p  
            for i in 1:Cp[p1].nc, j in 1:Cp[p2].nc 
                Avector=Cp[p1].Avector[i]+Cp[p2].Avector[j]   
                pos=findfirst(x->x==Avector,Cp[p].Avector)
                if Cp[p].corresp[pos]>0
                    Ψ₁ = Cp[p1].W[:,i] 
                    Ψ₂ = Cp[p2].W[:,j]
                    assembly_G!(MeshInfo,MatParams,Cp[p],pos, Ψ₁, Ψ₂)
                end 
            end
        end
    end     
end

function assembly_G!(MeshInfo::MeshShell,MatParams::MaterialParamShell,Cp::Parametrisation, entry::Int64, Ψ₁::Vector{ComplexF64}, Ψ₂::Vector{ComplexF64})
    mult = 1.0
    neq = MeshInfo.NFREEDOF
    F = zeros(ComplexF64, MeshInfo.NDOF)
    Fe = zeros(ComplexF64, 48)
    BL1 = zeros(ComplexF64, 6, 48)
    BL2 = zeros(ComplexF64, 6, 48)
    A = zeros(ComplexF64, 6, 9)
    dU = zeros(ComplexF64, 9)
    U1 = zeros(ComplexF64, MeshInfo.NDOF)
    U2 = zeros(ComplexF64, MeshInfo.NDOF)
    nhalf = Int(length(Ψ₁) ÷ 2)
    U1[MeshInfo.FREEDOF] = @view Ψ₁[1:nhalf]
    U2[MeshInfo.FREEDOF] = @view Ψ₂[1:nhalf]

    PreCalculateParam = MeshInfo.PreCalculateParam
    @inbounds for i = 1 : MeshInfo.NEle 
        eleDof = @view MeshInfo.ElementDOF[i, :]
        # extract element displacements
        U1e = collect(@view U1[eleDof])
        U2e = collect(@view U2[eleDof])
        # element internal force
        Integrate_G!(MatParams, PreCalculateParam ,U1e, U2e, i, Fe, BL1, BL2, A, dU)
        # assemble element force
        F[eleDof] .+= Fe * mult
    end
    # generate array from 1 to MeshInfo.NFREEDOF
    @inbounds for j = 1 : MeshInfo.NFREEDOF
        Cp.rhs[neq+MeshInfo.FREEDOF_list[j],entry] -= F[MeshInfo.FREEDOF[j]]
    end
end

function Integrate_G!(MatParams::MaterialParamShell, PreCalculate::PreCalculateParam, 
    U1e::Vector{ComplexF64}, U2e::Vector{ComplexF64}, ElementNum, 
    Fe::Vector{ComplexF64}, BL1::Matrix{ComplexF64}, BL2::Matrix{ComplexF64}, 
    A::Matrix{ComplexF64}, dU::Vector{ComplexF64})
    
    fill!(Fe, 0.0)
    fill!(BL1, 0.0)
    fill!(BL2, 0.0)
    fill!(A, 0.0)
    fill!(dU, 0.0)
    
    ξ, η, ζ, w, GassPointsNum = GassPointsRammeShell("2*2*2")
    
    σ₁ = Vector{ComplexF64}(undef, size(PreCalculate.B0[ElementNum, 1], 1))
    σ₂ = similar(σ₁)
    σ₁₂ = similar(σ₁)
    σ₂₁ = similar(σ₁)
    kau_nl = zeros(ComplexF64, 4, 48)
    kua_nl = zeros(ComplexF64, 48, 4)
    
    for i = 1:GassPointsNum
        StrainMatrixNonlinearQ8RammeShell!(BL1, A, dU, PreCalculate.G[ElementNum, i], U1e)
        StrainMatrixNonlinearQ8RammeShell!(BL2, A, dU, PreCalculate.G[ElementNum, i], U2e)
        
        B0 = PreCalculate.B0[ElementNum, i]
        D_i = PreCalculate.D[ElementNum, i]
        
        σ₁ = D_i * (B0 * U1e)
        σ₂ = D_i * (B0 * U2e)
        σ₁₂ = 0.5 * D_i * (BL1 * U2e)
        σ₂₁ = 0.5 * D_i * (BL2 * U1e)
        kau_nl += PreCalculate.Bα[ElementNum, i]'*PreCalculate.D[ElementNum, i]*BL1*PreCalculate.det_a[ElementNum, i]
        kua_nl += BL1'*PreCalculate.D[ElementNum, i]*PreCalculate.Bα[ElementNum, i]*PreCalculate.det_a[ElementNum, i]
        Fe .+= 0.5 * (B0'*σ₁₂ + B0'*σ₂₁ + BL1'*σ₂ + BL2'*σ₁) * PreCalculate.det_a[ElementNum, i]
    end
    Fe = Fe - 0.5 * PreCalculate.K0a[ElementNum, 1] *  PreCalculate.invKaa[ElementNum, 1] * kau_nl * U2e - kua_nl * PreCalculate.invKaa[ElementNum, 1] * PreCalculate.Ka0[ElementNum, 1] * U2e
end

function fillrhsH!(MeshInfo::MeshShell,MatParams::MaterialParamShell,Cp::Vector{Parametrisation},p::Int64)
    for p1 in 1:p-2, p2 in 1:p-2, p3 in 1:p-2
        if (p1+p2+p3)==p  
          for i in 1:Cp[p1].nc, j in 1:Cp[p2].nc, k in 1:Cp[p3].nc  
            Avector=Cp[p1].Avector[i]+Cp[p2].Avector[j]+Cp[p3].Avector[k]   
            pos=findfirst(x->x==Avector,Cp[p].Avector)
            if Cp[p].corresp[pos]>0
              Ψ₁ = Cp[p1].W[:,i]
              Ψ₂ = Cp[p2].W[:,j]
              Ψ₃ = Cp[p3].W[:,k]
              assembly_H!(MeshInfo,MatParams,Cp[p],pos, Ψ₁, Ψ₂, Ψ₃)
            end
          end
        end
    end     
end

function assembly_H!(MeshInfo::MeshShell,MatParams::MaterialParamShell,Cp::Parametrisation, entry::Int64, Ψ₁::Vector{ComplexF64}, Ψ₂::Vector{ComplexF64}, Ψ₃::Vector{ComplexF64})
    mult = 1.0
    neq = MeshInfo.NFREEDOF
    F = zeros(ComplexF64, MeshInfo.NDOF)
    U1 = zeros(ComplexF64, MeshInfo.NDOF)
    U2 = zeros(ComplexF64, MeshInfo.NDOF)
    U3 = zeros(ComplexF64, MeshInfo.NDOF)
    Fe = zeros(ComplexF64, 48)
    BL1 = zeros(ComplexF64, 6, 48)
    BL2 = zeros(ComplexF64, 6, 48)
    BL3 = zeros(ComplexF64, 6, 48)
    A = zeros(ComplexF64, 6, 9)
    dU = zeros(ComplexF64, 9)

    nhalf = Int(length(Ψ₁)/2)

    U1[MeshInfo.FREEDOF] = @view Ψ₁[1:nhalf]
    U2[MeshInfo.FREEDOF] = @view Ψ₂[1:nhalf]   
    U3[MeshInfo.FREEDOF] = @view Ψ₃[1:nhalf]
    PreCalculate = MeshInfo.PreCalculateParam
    @inbounds for i = 1 : MeshInfo.NEle
        eleDof = MeshInfo.ElementDOF[i, :]
        # extract element displacements
        U1e = collect(@view U1[eleDof])
        U2e = collect(@view U2[eleDof])
        U3e = collect(@view U3[eleDof])
        # element internal force
        Integrate_H!(MatParams, PreCalculate, U1e, U2e, U3e, i, Fe, BL1, BL2, BL3, A, dU)
        # assemble element force
        F[eleDof] .+= Fe * mult
    end 
    @inbounds for j = 1 : MeshInfo.NFREEDOF
        Cp.rhs[neq+MeshInfo.FREEDOF_list[j],entry] -= F[MeshInfo.FREEDOF[j]]
    end
end

function Integrate_H!(MatParams::MaterialParamShell, PreCalculate::PreCalculateParam, 
    U1e::Vector{ComplexF64}, U2e::Vector{ComplexF64}, U3e::Vector{ComplexF64}, ElementNum, 
    Fe::Vector{ComplexF64}, BL1::Matrix{ComplexF64}, BL2::Matrix{ComplexF64}, BL3::Matrix{ComplexF64}, 
    A::Matrix{ComplexF64}, dU::Vector{ComplexF64})
    
    fill!(Fe, 0.0)
    fill!(BL1, 0.0)
    fill!(BL2, 0.0)
    fill!(BL3, 0.0)
    fill!(A, 0.0)
    fill!(dU, 0.0)
    
    ξ, η, ζ, w, GassPointsNum = GassPointsRammeShell("2*2*2")
    kau_nl = zeros(ComplexF64, 4, 48)
    kua_nl = zeros(ComplexF64, 48, 4)
    
    for i = 1:GassPointsNum
        # 
        StrainMatrixNonlinearQ8RammeShell!(BL1, A, dU, PreCalculate.G[ElementNum, i], U1e)
        StrainMatrixNonlinearQ8RammeShell!(BL2, A, dU, PreCalculate.G[ElementNum, i], U2e)
        StrainMatrixNonlinearQ8RammeShell!(BL3, A, dU, PreCalculate.G[ElementNum, i], U3e)
        
        D_i = PreCalculate.D[ElementNum, i]
        
        # 
        σ₁₂ = 0.5 * D_i * (BL1 * U2e)
        σ₂₁ = 0.5 * D_i * (BL2 * U1e)
        σ₁₃ = 0.5 * D_i * (BL1 * U3e)
        σ₃₁ = 0.5 * D_i * (BL3 * U1e)
        σ₂₃ = 0.5 * D_i * (BL2 * U3e)
        σ₃₂ = 0.5 * D_i * (BL3 * U2e)
        kau_nl += PreCalculate.Bα[ElementNum, i]'*PreCalculate.D[ElementNum, i]*BL2*PreCalculate.det_a[ElementNum, i]
        kua_nl += BL1'*PreCalculate.D[ElementNum, i]*PreCalculate.Bα[ElementNum, i]*PreCalculate.det_a[ElementNum, i]
        Fe .+= (1/6) * (BL1'*σ₂₃ + BL1'*σ₃₂ + BL2'*σ₁₃ + BL2'*σ₃₁ + BL3'*σ₁₂ + BL3'*σ₂₁) * PreCalculate.det_a[ElementNum, i]
    end
    Fe = Fe - 0.5*kua_nl*PreCalculate.invKaa[ElementNum, 1]*kau_nl*U3e
end
