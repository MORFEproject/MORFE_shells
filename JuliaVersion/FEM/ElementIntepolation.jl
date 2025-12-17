function ShapeFunctionQ8Linear(s::Float64, t::Float64)
    N = zeros(Float64, 8)
    N[1] = 0.25 * (1 - s) * (1 - t) * (-s - t - 1)
    N[2] = 0.25 * (1 + s) * (1 - t) * (s - t - 1)
    N[3] = 0.25 * (1 + s) * (1 + t) * (s + t - 1)
    N[4] = 0.25 * (1 - s) * (1 + t) * (-s + t - 1)

    N[5] = 0.5 * (1 + s) * (1 - t) * (1 - s)
    N[6] = 0.5 * (1 + s) * (1 + t) * (1 - t)
    N[7] = 0.5 * (1 + s) * (1 + t) * (1 - s)
    N[8] = 0.5 * (1 - s) * (1 + t) * (1 - t)
    return N
end

function ShapeFunctionQ8Matrix(N::Array{Float64}, θ₃::Float64, h::Vector{Float64})
    Nmat1 = SubMatrix(N[1], θ₃, h[1])
    Nmat2 = SubMatrix(N[2], θ₃, h[2])
    Nmat3 = SubMatrix(N[3], θ₃, h[3])
    Nmat4 = SubMatrix(N[4], θ₃, h[4])
    Nmat5 = SubMatrix(N[5], θ₃, h[5])
    Nmat6 = SubMatrix(N[6], θ₃, h[6])
    Nmat7 = SubMatrix(N[7], θ₃, h[7])
    Nmat8 = SubMatrix(N[8], θ₃, h[8])
    Nmat = [Nmat1 Nmat2 Nmat3 Nmat4 Nmat5 Nmat6 Nmat7 Nmat8]
    
    return Nmat
end

function SubMatrix(N::Float64, θ₃::Float64, h::Float64)
    Nmat = zeros(Float64, 3, 6)
    Nmat[1,1] = N
    Nmat[2,2] = N
    Nmat[3,3] = N
    Nmat[1,4] = 0.5 * h * θ₃ * N
    Nmat[2,5] = 0.5 * h * θ₃ * N
    Nmat[3,6] = 0.5 * h * θ₃ * N
    return Nmat
end

function ShapeFunctionQ8LinearDerivative(s::Float64, t::Float64)
    dN = zeros(Float64, 8, 3)
    dN[1,1] = - (s/4 - 1/4)*(t - 1) - ((t - 1)*(s + t + 1))/4
    dN[1,2] = - (s/4 - 1/4)*(t - 1) - (s/4 - 1/4)*(s + t + 1)
    dN[1,3] = 0.0

    dN[2,1] = ((t - 1)*(t - s + 1))/4 - (s/4 + 1/4)*(t - 1)
    dN[2,2] = (s/4 + 1/4)*(t - s + 1) + (s/4 + 1/4)*(t - 1)
    dN[2,3] = 0.0

    dN[3,1] = (s/4 + 1/4)*(t + 1) + ((t + 1)*(s + t - 1))/4
    dN[3,2] = (s/4 + 1/4)*(t + 1) + (s/4 + 1/4)*(s + t - 1)
    dN[3,3] = 0.0  

    dN[4,1] = (s/4 - 1/4)*(t + 1) + ((t + 1)*(s - t + 1))/4
    dN[4,2] = (s/4 - 1/4)*(s - t + 1) - (s/4 - 1/4)*(t + 1)
    dN[4,3] = 0.0 

    dN[5,1] = (t/2 - 1/2)*(s - 1) + (t/2 - 1/2)*(s + 1)
    dN[5,2] = ((s - 1)*(s + 1))/2
    dN[5,3] = 0.0

    dN[6,1] = -((t - 1)*(t + 1))/2
    dN[6,2] = - (s/2 + 1/2)*(t - 1) - (s/2 + 1/2)*(t + 1)
    dN[6,3] = 0.0     

    dN[7,1] = - (t/2 + 1/2)*(s - 1) - (t/2 + 1/2)*(s + 1)
    dN[7,2] = -((s - 1)*(s + 1))/2
    dN[7,3] = 0.0     

    dN[8,1] = ((t - 1)*(t + 1))/2
    dN[8,2] = (s/2 - 1/2)*(t - 1) + (s/2 - 1/2)*(t + 1)
    dN[8,3] = 0.0  
    return dN   
end

function LocalBase(dN::Matrix{Float64}, N::Vector{Float64}, θ₃::Float64, X::Vector{Float64}, Y::Vector{Float64}, Z::Vector{Float64}, VectorDirection::Vector{Vector{Float64}}, h::Vector{Float64})
    a = zeros(Float64, 3, 3)
    for k = 1 : 8
        a[1,1] += dN[k,1]*X[k] + 0.5*θ₃*h[k]*dN[k,1]*VectorDirection[k][1]
        a[1,2] += dN[k,1]*Y[k] + 0.5*θ₃*h[k]*dN[k,1]*VectorDirection[k][2]
        a[1,3] += dN[k,1]*Z[k] + 0.5*θ₃*h[k]*dN[k,1]*VectorDirection[k][3]
        
        a[2,1] += dN[k,2]*X[k] + 0.5*θ₃*h[k]*dN[k,2]*VectorDirection[k][1]
        a[2,2] += dN[k,2]*Y[k] + 0.5*θ₃*h[k]*dN[k,2]*VectorDirection[k][2]
        a[2,3] += dN[k,2]*Z[k] + 0.5*θ₃*h[k]*dN[k,2]*VectorDirection[k][3]
        
        a[3,1] += 0.5*h[k]*N[k]*VectorDirection[k][1]
        a[3,2] += 0.5*h[k]*N[k]*VectorDirection[k][2]
        a[3,3] += 0.5*h[k]*N[k]*VectorDirection[k][3] 
    end
    return a
end

function ElasticMatrixRammeShell(E::Float64, ν::Float64, a::Array{Float64,2})
    aa = a*a'
    aa_inv = inv(aa)
    YG = E
    NU = ν
    LAMDA = YG*NU/((1.0 + NU)*(1.0-2.0*NU))
    MU = YG/(2.0*(1.0 + NU))
    LAMDA2MU = LAMDA+2.0*MU
    D=zeros(Float64, 6, 6)

    D[1,1]=LAMDA2MU*aa_inv[1,1]*aa_inv[1,1]
    D[1,2]=LAMDA*aa_inv[1,1]*aa_inv[2,2]+2.0*MU*aa_inv[1,2]*aa_inv[1,2]
    D[1,3]=LAMDA*aa_inv[1,1]*aa_inv[3,3]+2.0*MU*aa_inv[1,3]*aa_inv[1,3]
    D[1,4]=LAMDA2MU*aa_inv[1,1]*aa_inv[1,2]
    D[1,5]=LAMDA*aa_inv[1,1]*aa_inv[2,3]+2.0*MU*aa_inv[1,2]*aa_inv[1,3]
    D[1,6]=LAMDA2MU*aa_inv[1,1]*aa_inv[1,3]
    
    D[2,1]=D[1,2]
    D[2,2]=LAMDA2MU*aa_inv[2,2]*aa_inv[2,2]
    D[2,3]=LAMDA*aa_inv[2,2]*aa_inv[3,3]+2.0*MU*aa_inv[2,3]*aa_inv[2,3]
    D[2,4]=LAMDA2MU*aa_inv[2,2]*aa_inv[2,1]
    D[2,5]=LAMDA2MU*aa_inv[2,2]*aa_inv[2,3]
    D[2,6]=LAMDA*aa_inv[2,2]*aa_inv[3,1]+2.0*MU*aa_inv[2,3]*aa_inv[2,1]
    
    D[3,1]=D[1,3]
    D[3,2]=D[2,3]
    D[3,3]=LAMDA2MU*aa_inv[3,3]*aa_inv[3,3]
    D[3,4]=LAMDA*aa_inv[3,3]*aa_inv[1,2]+2.0*MU*aa_inv[3,1]*aa_inv[3,2]
    D[3,5]=LAMDA2MU*aa_inv[3,3]*aa_inv[3,2]
    D[3,6]=LAMDA2MU*aa_inv[3,3]*aa_inv[3,1]
    
    D[4,1]=D[1,4]
    D[4,2]=D[2,4]
    D[4,3]=D[3,4]
    D[4,4]=(LAMDA+MU)*aa_inv[1,2]*aa_inv[1,2]+MU*aa_inv[1,1]*aa_inv[2,2]
    D[4,5]=(LAMDA+MU)*aa_inv[1,2]*aa_inv[2,3]+MU*aa_inv[1,3]*aa_inv[2,2]
    D[4,6]=(LAMDA+MU)*aa_inv[1,2]*aa_inv[3,1]+MU*aa_inv[1,1]*aa_inv[2,3]
    
    D[5,1]=D[1,5]
    D[5,2]=D[2,5]
    D[5,3]=D[3,5]
    D[5,4]=D[4,5]
    D[5,5]=(LAMDA+MU)*aa_inv[2,3]*aa_inv[2,3]+MU*aa_inv[2,2]*aa_inv[3,3]
    D[5,6]=(LAMDA+MU)*aa_inv[2,3]*aa_inv[3,1]+MU*aa_inv[2,1]*aa_inv[3,3]
    
    D[6,1]=D[1,6]
    D[6,2]=D[2,6]
    D[6,3]=D[3,6]
    D[6,4]=D[4,6]
    D[6,5]=D[5,6]
    D[6,6]=(LAMDA+MU)*aa_inv[3,1]*aa_inv[3,1]+MU*aa_inv[3,3]*aa_inv[1,1]
    return D           
end

function CovarianteMatrix(a::Array{Float64,2})
    R=zeros(Float64, 6, 9)
    R[1,1:3]=a[1,1:3]
    R[2,4:6]=a[2,1:3]
    R[3,7:9]=a[3,1:3]
    
    R[4,1:3]=a[2,1:3]
    R[4,4:6]=a[1,1:3]
    
    R[5,4:6]=a[3,1:3]
    R[5,7:9]=a[2,1:3]
    
    R[6,1:3]=a[3,1:3]
    R[6,7:9]=a[1,1:3]
    return R
end

function GradientDeformation(θ₃::Float64, N::Vector{Float64}, dN::Matrix{Float64}, h::Vector{Float64})
    G=zeros(Float64, 9, 6*8)
    for  K = 1:8
        indice = 6*(K-1)
        TET3 = θ₃*0.5*h[K]
        G[1,1+indice]=dN[K,1]
        G[1,4+indice]=TET3*dN[K,1]
        G[2,2+indice]=G[1,1+indice]
        G[2,5+indice]=G[1,4+indice]
        G[3,3+indice]=G[1,1+indice]
        G[3,6+indice]=G[1,4+indice]
        G[4,1+indice]=dN[K,2]
        G[4,4+indice]=TET3*dN[K,2]
        G[5,2+indice]=G[4,1+indice]
        G[5,5+indice]=G[4,4+indice]
        G[6,3+indice]=G[4,1+indice]
        G[6,6+indice]=G[4,4+indice]
        G[7,4+indice]=0.5*h[K]*N[K]
        G[8,5+indice]=G[7,4+indice]
        G[9,6+indice]=G[7,4+indice]    
    end  
    return G
end

function StrainMatrixAlpha(θ₁::Float64, θ₂::Float64, θ₃::Float64, h::Vector{Float64}, N::Vector{Float64})
    epai=0;
    for k=1:8
        epai=epai+h[k]*N[k]
    end
    
    Bα=zeros(Float64, 6, 4)
    
    Bα[3,1]=epai/2*θ₃
    Bα[3,2]=epai/2*θ₃*θ₁
    Bα[3,3]=epai/2*θ₃*θ₂
    Bα[3,4]=epai/2*θ₃*θ₁*θ₂
    return Bα
end

function StrainMatrixNonlinearQ8RammeShell!(B::Matrix{ComplexF64}, A::Matrix{ComplexF64}, dU::Vector{ComplexF64}, G, U::AbstractVector{ComplexF64})
    mul!(dU, G, U)
    
    fill!(A, 0.0)
    @inbounds begin
        A[1,1:3] .= @view dU[1:3]
        A[2,4:6] .= @view dU[4:6]
        A[3,7:9] .= @view dU[7:9]
        A[4,1:3] .= @view dU[4:6]
        A[4,4:6] .= @view dU[1:3]
        A[5,4:6] .= @view dU[7:9]
        A[5,7:9] .= @view dU[4:6]
        A[6,1:3] .= @view dU[7:9]
        A[6,7:9] .= @view dU[1:3]
    end
    mul!(B, A, G)
    return B
end

function StrainMatrixNonlinearQ8RammeShell!(B::Matrix{Float64}, A::Matrix{Float64}, dU::Vector{Float64}, G, U::AbstractVector{Float64})
    mul!(dU, G, U)
    
    fill!(A, 0.0)
    @inbounds begin
        A[1,1:3] .= @view dU[1:3]
        A[2,4:6] .= @view dU[4:6]
        A[3,7:9] .= @view dU[7:9]
        A[4,1:3] .= @view dU[4:6]
        A[4,4:6] .= @view dU[1:3]
        A[5,4:6] .= @view dU[7:9]
        A[5,7:9] .= @view dU[4:6]
        A[6,1:3] .= @view dU[7:9]
        A[6,7:9] .= @view dU[1:3]
    end
    
    mul!(B, A, G)
    return B
end

function StressMatrixQ8Shell!(M::Matrix{Float64}, σ::Vector{Float64})
    fill!(M, 0.0)
    
    @inbounds begin
        for i = 1:3
            M[i, i] = σ[1]
            M[i+3, i+3] = σ[2]
            M[i+6, i+6] = σ[3]
        end
        
        for i = 1:3
            M[i, i+3] = σ[4]
            M[i+3, i] = σ[4]
            
            M[i, i+6] = σ[6]
            M[i+6, i] = σ[6]
            
            M[i+3, i+6] = σ[5]
            M[i+6, i+3] = σ[5]
        end
    end
end

# gass points for RammeShell elements
function GassPointsRammeShell(Type::String)
    if Type == "2*2*2"
        p = 1/sqrt(3);
        GassVector = [-p -p p;
                        p -p p;
                        p p p;
                        -p p p;
                        -p -p -p;
                        p -p -p;
                        p p -p;
                        -p p -p]
        ξ = GassVector[:,1]
        η = GassVector[:,2]
        ζ = GassVector[:,3]
        w = [1,1,1,1,1,1,1,1]
        GassPointsNum = 8

    elseif Type == "2*2*1"
        p = 1/sqrt(3)
        ξ = [-1.0, 1.0, -1.0, 1.0] * p
        η = [-1.0, -1.0, 1.0, 1.0] * p
        ζ = [0.0, 0.0, 0.0, 0.0]
        w = [1.0, 1.0, 1.0, 1.0]
        GassPointsNum = 4

    elseif Type == "2*2*3"
        p1 = 1/sqrt(3)
        p2 = sqrt(3/5)
        w1 = 8/9
        w2 = 5/9
        ξ = [-1.0, -1.0, -1.0,  1.0,  1.0,  1.0, -1.0, -1.0, -1.0,  1.0,  1.0,  1.0] * p1
        η = [-1.0, -1.0, -1.0, -1.0, -1.0, -1.0,  1.0,  1.0,  1.0,  1.0,  1.0,  1.0] * p1
        ζ = [-1.0,  0.0,  1.0, -1.0,  0.0,  1.0, -1.0,  0.0,  1.0, -1.0,  0.0,  1.0] * p2
        w =  [w2, w1, w2, w2, w1, w2, w2, w1, w2, w2, w1, w2]
        GassPointsNum = 12

    elseif Type == "3*3*2"
        p1 = sqrt(3/5)
        p2 = 1/sqrt(3)
        w1 = 8/9
        w2 = 5/9
        ξ = [-1.0, -1.0, 0.0, 0.0, 1.0, 1.0, -1.0, -1.0, 0.0, 0.0, 1.0, 1.0, -1.0, -1.0, 0.0, 0.0, 1.0, 1.0] * p1
        η = [-1.0, -1.0, -1.0, -1.0, -1.0, -1.0,  0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0,  1.0,  1.0,  1.0,  1.0,  1.0] * p1
        ζ = [-1.0,  1.0, -1.0,  1.0, -1.0,  1.0, -1.0,  1.0, -1.0,  1.0, -1.0,  1.0, -1.0,  1.0, -1.0,  1.0, -1.0,  1.0] * p2
        wp1 = [w2;w2;w1;w1;w2;w2;w2;w2;w1;w1;w2;w2;w2;w2;w1;w1;w2;w2]
        wp2 = [w2;w2;w2;w2;w2;w2;w1;w1;w1;w1;w1;w1;w2;w2;w2;w2;w2;w2]
        w = wp1 .* wp2
        GassPointsNum = 18

    elseif Type == "3*3*3"
        p = 1/sqrt(3);
        GpLoc = [-p -p p;
                    p -p p;
                    p p p;
                    -p p p;
                    -p -p -p;
                    p -p -p;
                    p p -p;
                    -p p -p; # 1-8
                    p p 0;
                    p -p 0;
                    -p p 0;
                    -p -p 0; # 9-12
                    p 0 p;
                    p 0 -p;
                    -p 0 p;
                    -p 0 -p; # 13-16
                    0 p p;
                    0 p -p;
                    0 -p p;
                    0 -p -p; # 17-20
                    0 0 p;
                    0 0 -p; # 21-22
                    0 p 0;
                    0 -p 0; # 23-24
                    p 0 0;
                    -p 0 0; # 25-26
                    0 0 0] # 27
            ξ = GpLoc[:,1]
            η = GpLoc[:,2]
            ζ = GpLoc[:,3]
            w = [(5/9)^3, (5/9)^3, (5/9)^3, (5/9)^3, (5/9)^3, (5/9)^3, (5/9)^3, (5/9)^3, 
                (5/9)^2*(8/9), (5/9)^2*(8/9), (5/9)^2*(8/9), (5/9)^2*(8/9), 
                (5/9)^2*(8/9), (5/9)^2*(8/9), (5/9)^2*(8/9), (5/9)^2*(8/9), 
                (5/9)^2*(8/9), (5/9)^2*(8/9), (5/9)^2*(8/9), (5/9)^2*(8/9), 
                (5/9)*(8/9)^2, (5/9)*(8/9)^2, (5/9)*(8/9)^2, (5/9)*(8/9)^2, (5/9)*(8/9)^2, (5/9)*(8/9)^2, 
                (8/9)^3]
            GassPointsNum = 27;
    end
    return ξ, η, ζ, w, GassPointsNum
end