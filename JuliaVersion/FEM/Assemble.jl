function AssembleKMRammeShellQ8(MeshInfo::MeshShell, MatParams::MaterialParamShell)
    K = spzeros(Float64, MeshInfo.NDOF, MeshInfo.NDOF)
    M = spzeros(Float64, MeshInfo.NDOF, MeshInfo.NDOF)
    for i = 1 : MeshInfo.NEle
        Kₑ = MeshInfo.PreCalculateParam.K0[i] - MeshInfo.PreCalculateParam.K0a[i] * MeshInfo.PreCalculateParam.invKaa[i] * MeshInfo.PreCalculateParam.Ka0[i]
        Mₑ = MeshInfo.PreCalculateParam.M[i]
        eleDof = MeshInfo.ElementDOF[i, :]
        K[eleDof, eleDof] += Kₑ
        M[eleDof, eleDof] += Mₑ
    end
    return K, M
end
# total lagrange pre-calculation
function PreCalculate!(MeshInfo::MeshShell, MatParams::MaterialParamShell)
    ξₖ, ηₖ, ζₖ, wₖ, GassPointsNumₖ = GassPointsRammeShell("2*2*2") # 2*2*2
    ξₘ, ηₘ, ζₘ, wₘ, GassPointsNumₘ = GassPointsRammeShell("3*3*2") # 3*3*2
    det_a = zeros(Float64, MeshInfo.NEle, GassPointsNumₖ)
    D = Array{Any}(undef, MeshInfo.NEle, GassPointsNumₖ) # elastic matrix
    R = Array{Any}(undef, MeshInfo.NEle, GassPointsNumₖ) # covariant tensor matrix for each element integration point
    B0 = Array{Any}(undef, MeshInfo.NEle, GassPointsNumₖ) # linear strain matrix for each element integration point
    G = Array{Any}(undef, MeshInfo.NEle, GassPointsNumₖ)
    Bα = Array{Any}(undef, MeshInfo.NEle, GassPointsNumₖ) # enhanced strain matrix for each element integration point

    invKaa = Array{Any}(undef, MeshInfo.NEle) # inverse enhanced strain stiffness matrix for each element
    K0 = Array{Any}(undef, MeshInfo.NEle) # linear stiffness matrix for each element
    K0a = Array{Any}(undef, MeshInfo.NEle) # inverse diagonal matrix -1 for each element
    Ka0 = Array{Any}(undef, MeshInfo.NEle) # inverse diagonal matrix -2 for each element
    M = Array{Any}(undef, MeshInfo.NEle) # mass matrix for each element

    E = MatParams.E
    ν = MatParams.ν
    ρ = MatParams.ρ
    

    for i = 1 : MeshInfo.NEle
        ElementConnectivity = MeshInfo.ElementNode[i, :]
        X = MeshInfo.NodeCoord[ElementConnectivity, 1]
        Y = MeshInfo.NodeCoord[ElementConnectivity, 2]
        Z = MeshInfo.NodeCoord[ElementConnectivity, 3]
        e1 = MeshInfo.ElementNormal[i, 1]
        e2 = MeshInfo.ElementNormal[i, 2]
        e3 = MeshInfo.ElementNormal[i, 3]
        e4 = MeshInfo.ElementNormal[i, 4]
        e5 = MeshInfo.ElementNormal[i, 5]
        e6 = MeshInfo.ElementNormal[i, 6]
        e7 = MeshInfo.ElementNormal[i, 7]
        e8 = MeshInfo.ElementNormal[i, 8]
        VectorDirection = [e1, e2, e3, e4, e5, e6, e7, e8]
        t = MatParams.t[ElementConnectivity]
        kaa = zeros(Float64, 4, 4)
        k0 = zeros(Float64, 48, 48)
        k0a = zeros(Float64, 48, 4)
        ka0 = zeros(Float64, 4, 48)
        m = zeros(Float64, 48, 48)

        for j = 1 : GassPointsNumₖ
            θ₁ = ξₖ[j]
            θ₂ = ηₖ[j]
            θ₃ = ζₖ[j]
            N = ShapeFunctionQ8Linear(θ₁, θ₂)
            dN = ShapeFunctionQ8LinearDerivative(θ₁, θ₂)
            a = LocalBase(dN, N, θ₃, X, Y, Z, VectorDirection, t)
            det_a[i, j] = det(a) * wₖ[j]
            D[i, j] = ElasticMatrixRammeShell(E, ν, a)
            R[i, j] = CovarianteMatrix(a)
            G[i, j] = GradientDeformation(θ₃, N, dN, t)
            Bα[i, j] = StrainMatrixAlpha(θ₁, θ₂, θ₃, t, N)
            B0[i, j] = R[i, j] * G[i, j]
            kaa += Bα[i, j]' * D[i, j] * Bα[i, j] * det_a[i, j]
            k0 += B0[i, j]' * D[i, j] * B0[i, j] * det_a[i, j]
            k0a += B0[i, j]' * D[i, j] * Bα[i, j] * det_a[i, j]
            ka0 += Bα[i, j]' * D[i, j] * B0[i, j] * det_a[i, j]
        end

        for k = 1 : GassPointsNumₘ
            θ₁ = ξₘ[k]
            θ₂ = ηₘ[k]
            θ₃ = ζₘ[k]
            N = ShapeFunctionQ8Linear(θ₁, θ₂)
            Nmat = ShapeFunctionQ8Matrix(N, θ₃, t)
            dN = ShapeFunctionQ8LinearDerivative(θ₁, θ₂)
            a = LocalBase(dN, N, θ₃, X, Y, Z, VectorDirection, t)
            m += Nmat' * ρ * Nmat * det(a) * wₘ[k]
        end
        invKaa[i] = inv(kaa)
        K0[i] = k0
        K0a[i] = k0a
        Ka0[i] = ka0
        M[i] = m
    end
    MeshInfo.PreCalculateParam.invKaa = invKaa
    MeshInfo.PreCalculateParam.K0 = K0
    MeshInfo.PreCalculateParam.K0a = K0a
    MeshInfo.PreCalculateParam.Ka0 = Ka0
    MeshInfo.PreCalculateParam.M = M
    MeshInfo.PreCalculateParam.det_a = det_a

    MeshInfo.PreCalculateParam.D = D
    MeshInfo.PreCalculateParam.R = R
    MeshInfo.PreCalculateParam.B0 = B0
    MeshInfo.PreCalculateParam.G = G
    MeshInfo.PreCalculateParam.Bα = Bα

end

function ThicknessToGassPoint(s, t, he)
    N = ShapeFunctionQ4Linear(s, t)
    hp = N[1]*he[1] + N[2]*he[2] + N[3]*he[3] + N[4]*he[4]
    return hp
end


function StiffnessMatrixNonlinear!(KL::Matrix{Float64}, BNL::Matrix{Float64},
    A::Matrix{Float64}, dU::Vector{Float64}, temp1::Matrix{Float64}, temp2::Matrix{Float64},
    U::Vector{Float64}, MeshInfo::MeshShell, CurrentElement::Int)
    
    fill!(KL, 0.0)
    ξ, η, ζ, w, GassPointsNum = GassPointsRammeShell("2*2*2")
    
    @fastmath @inbounds for i = 1:GassPointsNum
        G = MeshInfo.PreCalculateParam.G[CurrentElement, i]
        B0 = MeshInfo.PreCalculateParam.B0[CurrentElement, i]
        StrainMatrixNonlinearQ8RammeShell!(BNL, A, dU, G, U)
        det_a = MeshInfo.PreCalculateParam.det_a[CurrentElement, i]
        D = MeshInfo.PreCalculateParam.D[CurrentElement, i]
        
        mul!(temp1, B0', D)
        mul!(temp2, temp1, BNL)
        @simd for j in eachindex(temp2)
            KL[j] += 0.5 * temp2[j] * det_a
        end
        
        mul!(temp1, BNL', D)
        mul!(temp2, temp1, B0)
        @simd for j in eachindex(temp2)
            KL[j] += temp2[j] * det_a
        end
        
        mul!(temp1, BNL', D)
        mul!(temp2, temp1, BNL)
        @simd for j in eachindex(temp2)
            KL[j] += 0.5 * temp2[j] * det_a
        end
    end
end

function StiffnessMatrixTangent!(KT::Matrix{Float64}, BNL::Matrix{Float64},
    A::Matrix{Float64}, dU::Vector{Float64}, temp1::Matrix{Float64}, temp2::Matrix{Float64},
    S_mat::Matrix{Float64}, σTemp::Vector{Float64}, B0_BNL::Matrix{Float64}, temp_strain::Matrix{Float64},
    temp1_sigam::Matrix{Float64},
    U::Vector{Float64}, MeshInfo::MeshShell, CurrentElement::Int)
    
    fill!(KT, 0.0)
    ξ, η, ζ, w, GassPointsNum = GassPointsRammeShell("2*2*2")
    
    @fastmath @inbounds for i = 1:GassPointsNum
        G = MeshInfo.PreCalculateParam.G[CurrentElement, i]
        B0 = MeshInfo.PreCalculateParam.B0[CurrentElement, i]
        StrainMatrixNonlinearQ8RammeShell!(BNL, A, dU, G, U)
        det_a = MeshInfo.PreCalculateParam.det_a[CurrentElement, i]
        D = MeshInfo.PreCalculateParam.D[CurrentElement, i]
        fill!(B0_BNL, 0.0)
        fill!(temp1, 0.0)
        @simd for j in 1:size(B0, 1)
            @simd for k in 1:size(B0, 2)
                B0_BNL[j, k] = B0[j, k] + 0.5 * BNL[j, k]
            end
        end
        
        fill!(temp_strain, 0.0)
        mul!(temp_strain, D, B0_BNL)      
        mul!(σTemp, temp_strain, U) 
        
        StressMatrixQ8Shell!(S_mat, σTemp)
        
        mul!(temp1, B0', D) # 48*6
        mul!(temp2, temp1, BNL) # 48*48
        @simd for j in eachindex(temp2)
            KT[j] += temp2[j] * det_a
        end
        
        mul!(temp1, BNL', D) # 48*6
        mul!(temp2, temp1, B0) # 48*48
        @simd for j in eachindex(temp2)
            KT[j] += temp2[j] * det_a
        end
        
        mul!(temp1, BNL', D) # 48*6
        mul!(temp2, temp1, BNL) # 48*48
        @simd for j in eachindex(temp2)
            KT[j] += temp2[j] * det_a
        end
        
        fill!(temp1_sigam, 0.0)
        mul!(temp1_sigam, G', S_mat) # 48*9
        mul!(temp2, temp1_sigam, G) # 48*48
        @simd for j in eachindex(temp2)
            KT[j] += temp2[j] * det_a
        end
    end
end

function mass_normalization!(ϕ,M,neig)
  # Normalize the eigenvectors in the mass matrix
  for i = 1:neig
    c = transpose(ϕ[:,i])*M*ϕ[:,i]
    for j = 1:M.m
      ϕ[j,i] /= sqrt(c)
    end
  end
  return nothing
end