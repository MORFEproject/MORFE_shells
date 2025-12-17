"""
populate_with_complex_representation!(Cp::Vector{Parametrisation}, info::SSMParam)

This function populates the 'real' coefficient matrices (Wr, fr) directly with 
  the coefficients from the 'complex' representation (W, f). It serves as an alternative to realification!, as it does not perform a coordinate basis transformation. Instead, it effectively retains the complex coordinate basis for all subsequent calculations.

This approach assumes that downstream functions are aware that these coefficients 
  correspond to the complex variables (z_k, z̄_k), and not the real variables (x_k, y_k).
"""
function populate_with_complex_representation!(Cp::Vector{Parametrisation}, info::DPIMParam)
    for p in 1:info.max_order
        if size(Cp[p].Wr) != size(Cp[p].W) || size(Cp[p].fr) != size(Cp[p].f)
            error("dimension mismatch: Wr/fr and W/f sizes are inconsistent at p=$p.")
        end
        Cp[p].Wr .= Cp[p].W
        Cp[p].fr .= Cp[p].f
    end
end

function realification!(Cp::Vector{Parametrisation},info::DPIMParam)

  for p in 1:info.max_order
    for i in 1:Cp[p].nc
      Avec=Cp[p].Avector[i][:]
      Ivec=zeros(Int64,p)
      counter=0
      for j in 1:info.nrom 
        Ivec[counter+1:counter+Avec[j]].=j
        counter+=Avec[j]
      end
      pos=1
      coeff=1.0+0.0im
      Avec.=0
      recursive_C2R!(Ivec,p,pos,i,Avec,coeff,Cp[p],info)
    end  
  end

end  

function recursive_C2R!(Ivec::Vector{Int64},p::Int64,pos::Int64,posinit::Int64,Avec::Vector{Int64},
                             coeff::ComplexF64,Cp::Parametrisation,info::DPIMParam)

  nzhalf=Int(info.nz/2)
  nzfhalf=Int(info.nzforce/2)                           
  if pos==(p+1)
    pos1=findfirst(x->x==Avec,Cp.Avector)
    Cp.Wr[:,pos1]+=coeff*Cp.W[:,posinit]
    Cp.fr[:,pos1]+=coeff*Cp.f[:,posinit]
  else
    Avec1=Avec[:]   # new vectors
    Avec2=Avec[:]
    iz=Ivec[pos]    # z var 

    if iz <= nzhalf
      coeff1=0.5*coeff  
      Avec1[iz]+=1
      coeff2=-0.5*im*coeff
      Avec2[iz+nzhalf]+=1
    elseif iz <= info.nz
      coeff1=0.5*coeff
      Avec1[iz-nzhalf]+=1
      coeff2=0.5*im*coeff
      Avec2[iz]+=1
    elseif iz <= info.nz+nzfhalf  # to have cos and sin the 1/2 must avoided
      coeff1=im*coeff
      Avec1[iz]+=1
      coeff2=coeff
      Avec2[iz+nzfhalf]+=1
    else    
      coeff1=-im*coeff
      Avec1[iz-nzfhalf]+=1
      coeff2=coeff
      Avec2[iz]+=1
    end     

    pos+=1
    recursive_C2R!(Ivec,p,pos,posinit,Avec1,coeff1,Cp,info)
    recursive_C2R!(Ivec,p,pos,posinit,Avec2,coeff2,Cp,info)
  end

end    
