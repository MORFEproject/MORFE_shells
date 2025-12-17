using Printf
# Homological equation solver for half state space equation
function homological_HALF!(Sol::Vector{ComplexF64}, Rhs::Vector{ComplexF64}, Mat::SparseMatrixCSC{ComplexF64},
                           Λ::Vector{ComplexF64}, Cp::Vector{Parametrisation}, p::Int64, Apos::Int64,
                           M::SparseMatrixCSC{Float64}, K::SparseMatrixCSC{Float64}, C::SparseMatrixCSC{Float64},
                           AY::Matrix{ComplexF64}, info::SSMParam, MatParams::MaterialParamShell,
                           logfile::IOStream)

    σ = dot(Cp[p].Avector[Apos], Λ[1:info.nrom])
    resonant_modes = zeros(Bool, info.nrom)
    # To use the CNF option, import the frequency combinations into the log file
    if info.style == "c"
        fill!(resonant_modes, false)
        for i = 1:info.nz
            λ = Λ[i]
              if (abs(imag(σ - λ)) / abs(imag(λ)) <= 1e-2)
                  resonant_modes[i] = true

                  # 1. judge resonance type
                  resonance_type_str = "Internal Resonance"
                  check_vector = copy(Cp[p].Avector[Apos])
                  check_vector[i] -= 1
                  if length(check_vector) == 2 * info.nm
                      simplified_coeffs = [check_vector[k] - check_vector[k + info.nm] for k in 1:info.nm]
                      if all(isapprox.(simplified_coeffs, 0; atol=1e-9))
                          resonance_type_str = "Trivial Resonance"
                      end
                  end

                  # 2. construct the string representation of sigma (using the new conjugate representation)
                  sigma_str_builder = IOBuffer()
                  first_term = true
                  for j in 1:length(Cp[p].Avector[Apos])
                      coeff = Cp[p].Avector[Apos][j]
                      if coeff == 0; continue; end

                      if !first_term
                          if coeff > 0; print(sigma_str_builder, " + "); else; print(sigma_str_builder, " - "); coeff = -coeff; end
                      elseif coeff < 0; print(sigma_str_builder, "-"); coeff = -coeff; end
                      
                      if coeff != 1; print(sigma_str_builder, coeff, "*"); end

                      if j <= info.nm
                          print(sigma_str_builder, "λ", j)
                      else
                          primary_index = j - info.nm
                          print(sigma_str_builder, "conj(λ", primary_index, ")")
                      end
                      first_term = false
                  end
                  sigma_str = String(take!(sigma_str_builder))
                  
                  lambda_i_str = i <= info.nm ? "λ$i" : "conj(λ$(i - info.nm))"

                  # 3. log the information
                  if resonance_type_str == "Internal Resonance"
                    @printf(logfile, "Detected: %s\n", "Internal Resonance")
                    @printf(logfile, "Polynomial Order: %d\n", p)
                    @printf(logfile, "Frequency Check: %s ≈ %s\n\n", sigma_str, lambda_i_str)
                  end

                  # log the information
                  println("--------------------------------------------------")
                  @printf("--> Detected: %s (Logged to file)\n", resonance_type_str)
                  @printf("    Frequency Check: %s ≈ λ%d\n", sigma_str, i)
                  @printf("    Relative Difference: %.4f\n", abs(imag(σ-λ))/abs(imag(λ)))
                  println("--------------------------------------------------")
              end
        end
    elseif info.style == "g"
        fill!(resonant_modes, true)
    end
    println(Cp[p].Avector[Apos],"  ",resonant_modes)
    fill!(Rhs,0.0)
    Rhs[1:info.nK]=Cp[p].rhs[info.nK+1:info.nA,Apos]-M*Cp[p].Wf[info.nK+1:info.nA,Apos]-(σ*M+C)*Cp[p].Wf[1:info.nK,Apos]
    fill!(Mat.nzval,0.0)
    Mat[1:info.nK,1:info.nK] = (σ^2+MatParams.α*σ)*M + (MatParams.β*σ+1)*K
    for d = 1:info.nz
    Mat[info.nK+d,info.nK+d] = 1.0       
    if resonant_modes[d]
    λ = Λ[d]
    Mat[info.nK+d,1:info.nK] = (σ-conj(λ))*AY[1:info.nK,d] 
    Mat[1:info.nK,info.nK+d] = (σ-conj(λ))*AY[1:info.nK,d] 
    if d<=info.nm && resonant_modes[d+info.nm]  
    Mat[info.nK+d,info.nK+d+info.nm] = 1.0 
    end  
    if d>info.nm && resonant_modes[d-info.nm]  
    Mat[info.nK+d,info.nK+d-info.nm] = 1.0 
    end  
    Rhs[info.nK+d]=-transpose(AY[1:info.nK,d])*Cp[p].Wf[1:info.nK,Apos]
    end  
    end  
    Sol.=Mat\Rhs
    Cp[p].W[1:info.nK,Apos]=Sol[1:info.nK]
    Cp[p].W[info.nK+1:info.nA,Apos]=σ*Sol[1:info.nK]+Cp[p].Wf[1:info.nK,Apos]
    Cp[p].f[1:info.nz,Apos]=Sol[info.nK+1:info.nK+info.nz]
    for d = 1:info.nz
    if resonant_modes[d]
    Cp[p].W[info.nK+1:info.nA,Apos]+=Cp[p].f[d,Apos]*Cp[1].W[1:info.nK,d]
    end  
    end
end

# Homological equation solver for full state space equation
function homological_FULL!(Sol::Vector{ComplexF64},Rhs::Vector{ComplexF64},Mat::SparseMatrixCSC{ComplexF64},
                           Λ::Vector{ComplexF64},Cp::Parametrisation,Apos::Int64,
                           M::SparseMatrixCSC{Float64},K::SparseMatrixCSC{Float64},C::SparseMatrixCSC{Float64},
                           AY::Matrix{ComplexF64},XTA::Matrix{ComplexF64},info::SSMParam)

                            
  σ=dot(Cp.Avector[Apos],Λ[1:info.nrom])
  resonant_modes = zeros(Bool,info.nrom)
  if info.style=="c"
    fill!(resonant_modes,false)
    for i = 1:info.nz
        λ_homo = Λ[i]  
        if (abs((σ-λ_homo))/abs((λ_homo))<=1e-3)
            resonant_modes[i] = true
        end
    end
  elseif (info.style=="g")
      fill!(resonant_modes,true)
  elseif (info.style=="r")
      fill!(resonant_modes,false)
      for i = 1:info.nz
          λ_homo = Λ[i] 
          if (abs(imag(σ-λ_homo))/abs(imag(λ_homo))<=1e-3)
              resonant_modes[i] = true
              if (i <= Int(info.nrom/2))
                  resonant_modes[i+Int(info.nrom/2)] = true
              else
                  resonant_modes[i-Int(info.nrom/2)] = true
              end
          end
      end
  end
  println(Cp.Avector[Apos],"  ",resonant_modes)

  fill!(Rhs,0.0)
  Rhs[1:info.nA]=Cp.rhs[:,Apos]
  Rhs[1:info.nK]-=M*Cp.Wf[1:info.nK,Apos]
  Rhs[info.nK+1:info.nA]-=M*Cp.Wf[info.nK+1:info.nA,Apos]

  fill!(Mat.nzval,0.0)
  Mat[1:info.nK,1:info.nK] = σ*M
  Mat[info.nK+1:info.nA,info.nK+1:info.nA] = σ*M+C
  Mat[1:info.nK,info.nK+1:info.nA] = -M
  Mat[info.nK+1:info.nA,1:info.nK] = K
  for j = 1:info.nz
    if resonant_modes[j]  # if resonant
      Mat[1:info.nA,info.nA+j]=AY[:,j]
      Mat[info.nA+j,1:info.nA]=XTA[j,:]
    else
      Mat[info.nA+j,info.nA+j]=1   
    end  
  end  
  Sol.=Mat\Rhs
  Cp.W[:,Apos]=Sol[1:info.nA]
  Cp.f[1:info.nz,Apos]=Sol[info.nA+1:info.nA+info.nz]

end

function fillWf!(Cp::Vector{Parametrisation},p::Int64,info::SSMParam)

  for p1 in 2:p-1, p2 in 2:p-1
    if (p1+p2)==p+1  
      for i in 1:Cp[p1].nc
        A1=Cp[p1].Avector[i][:]
        for j in 1:Cp[p2].nc 
          A2=Cp[p2].Avector[j][:]      
          for s in 1:info.nrom
            if A1[s]>0
              Avector=A1+A2
              Avector[s]-=1
              pos=findfirst(x->x==Avector,Cp[p].Avector)
              Cp[p].Wf[1:info.nK,pos]+=A1[s]*Cp[p1].W[1:info.nK,i]*Cp[p2].f[s,j]
              Cp[p].Wf[info.nK+1:info.nA,pos]+=A1[s]*Cp[p1].W[info.nK+1:info.nA,i]*Cp[p2].f[s,j]
            end  
          end    
        end 
      end
    end
  end

end 
    

function fillWfnonaut!(Cp::Vector{Parametrisation},p::Int64,Avector::Vector{Int64},ind_rhs::Int64,info::SSMParam)

  for r in info.nz+1:info.nrom
    if Avector[r]>0
      for s in 1:info.nz
        fs_r = Cp[1].f[s,r]
        if abs(fs_r)>10^(-8)    # only fills if the reduced dyn is nzero
          Av_W=Avector[:]
          Av_W[s]+=1
          Av_W[r]-=1
          ind_W=findfirst(x->x==Av_W,Cp[p].Avector)
          Cp[p].Wf[1:info.nK,ind_rhs]+=Cp[p].W[1:info.nK,ind_W]*fs_r*Av_W[s]
          Cp[p].Wf[info.nK+1:info.nA,ind_rhs]+=Cp[p].W[info.nK+1:info.nA,ind_W]*fs_r*Av_W[s]
        end
      end
    end  
  end

end

