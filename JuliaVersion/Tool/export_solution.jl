function count_terms_dyn(Cp::Vector{Parametrisation},info::DPIMParam)
    # Routine that evaluates how many terms are present
    howmany=0
    for p in 1:info.max_order   # for every order
        for i in 1:Cp[p].nc # for every alpha vector
           howmany=howmany+1
        end
    end
    return howmany
end

function store_dyn_and_map(Cp::Vector{Parametrisation},info::DPIMParam,howmany::Integer, mesh::MeshShell, MatParams::MaterialParamShell)
  mappings=zeros(howmany,mesh.NNode,6)
  mappings_vel=zeros(howmany,mesh.NNode,6)
  mappings_modal=zeros(howmany,info.neig)
  mappings_modal_vel=zeros(howmany,info.neig)
  Avector=zeros(howmany,Cp[1].nc)
  fdyn=zeros(howmany,info.nm*2)
  KFULL, MFULL = AssembleKMRammeShellQ8(mesh, MatParams)
  K = KFULL[mesh.FREEDOF, mesh.FREEDOF]
  M = MFULL[mesh.FREEDOF, mesh.FREEDOF]
  neig = maximum([info.Φ; info.neig]) 
  λ, ϕ = eigs(K,M,nev=neig,which=:SM)
  λ = real(λ)
  ϕ = real(ϕ)
  println(size(ϕ))
  ω = sqrt.(λ)
  println("ω = ", ω/2/pi)  

  PhiFull = zeros(mesh.NDOF, neig)
  PhiFull[mesh.FREEDOF, :] = ϕ
  
  for i in 1:neig
      PhiSelectFull = PhiFull[:, i]
      ncols = div(length(PhiSelectFull), 6)
      PhiSelectFull = reshape(PhiSelectFull, 6, ncols)
      if sum(PhiSelectFull[1:3, :]) < 0 || abs(minimum(PhiSelectFull[1:3, :])) > abs(maximum(PhiSelectFull[1:3, :]))
          ϕ[:, i] .= -ϕ[:, i]
      end
  end
  index=1
  for p in 1:info.max_order   # for every order
    for i in 1:Cp[p].nc # for every alpha vector
        Avector[index,:]=Cp[p].Avector[i]
        index=index+1
    end
  end

  index=1
  for p in 1:info.max_order   # for every order
    for i in 1:Cp[p].nc # for every alpha vector
        for j in 1:info.nm
          fdyn[index,j]=2*real(Cp[p].fr[j,i])
        end  
        for j in 1:info.nm
          fdyn[index,j+info.nm]=-2*imag(Cp[p].fr[j,i])
        end
        index=index+1   
    end
  end

  index=1
  Dof = FindDof(mesh.NDOF, mesh.FREEDOF)
  for p in 1:info.max_order   # for every order
    for i in 1:Cp[p].nc # for every alpha vector
        for inode=1:mesh.NNode # loop over the nodes
            for idof=1:6 # loop over the dof
              dof = Dof[6*(inode-1)+idof];
              if dof>0 
                mappings[index,inode,idof]=real(Cp[p].Wr[dof,i])
                mappings_vel[index,inode,idof]=real(Cp[p].Wr[info.nK+dof,i])
              end
            end
        end
        for imode=1:info.neig # loop over the modes
          mappings_modal[index,imode]=ϕ[:,imode]'*M*real.(Cp[p].Wr[1:info.nK,i])
          mappings_modal_vel[index,imode]=ϕ[:,imode]'*M*real.(Cp[p].Wr[info.nK+1:2*info.nK,i])
        end
        index=index+1
    end      
  end
  return  mappings,mappings_vel,mappings_modal,mappings_modal_vel,Avector,fdyn
end

function FindDof(sdof, freedof)
  freedof_list = 1:length(freedof)
  dof = fill(-1, sdof)
  dof[freedof] = freedof_list
  return dof
end



function write_rdyn(info::DPIMParam,Cp::Vector{Parametrisation})
    println("========================================")
    println("Begin write Reduced dynamic equations!")
    println("norm form = ", info.nrom )
    rdyn = ["" for i in 1:info.nz]
    for i = 1:info.nz
      rdyn[i] = "z"*string(i)*"' = "
    end
    
    for p in 1:info.max_order
      for c = 1:Cp[p].nc
        Avector=Cp[p].Avector[c]   
        monomial = ""
        for d = 1:info.nrom        
          if (Avector[d]!=0)
            monomial *= "*z"*string(d)*"^"*string(Avector[d])
          end
        end
        for j in 1:info.nm
          rcoeff=2*real(Cp[p].fr[j,c])
          icoeff=-2*imag(Cp[p].fr[j,c])
          if abs(rcoeff)>1e-20
            rdyn[j] *= " + "*string(rcoeff)*monomial
          end
          if abs(icoeff)>1e-20
            rdyn[j+info.nm] *= " + "*string(icoeff)*monomial
          end
        end
      end
    end
    for i = 1 : info.nz
        println(rdyn[i])
    end
    ofile = open("./Output/Equations.txt","w")
    for i = 1:info.nz
      write(ofile,rdyn[i]*";\n")
    end
    close(ofile)  
    println("Ending of writing equation!")
    println("========================================")
  end

"""
write_rdyn_from_complex(info::DPIMParam, Cp::Vector{Parametrisation})
"""
function write_rdyn_from_complex(info::DPIMParam, Cp::Vector{Parametrisation})
    println("========================================")
    println("Begin write Reduced dynamic equations!")
    println("norm form = ", info.nrom )
    nm = info.nm

    rdyn_z = ["" for i in 1:nm]
    for i in 1:nm
        rdyn_z[i] = "z" * string(i) * "' = "
    end

    for p in 1:info.max_order
        for c in 1:Cp[p].nc
            Avector = Cp[p].Avector[c]
            monomial_str = ""
            for d in 1:info.nrom 
                if Avector[d] != 0
                    var_name = ""
                    if d <= nm
                        var_name = "z" * string(d)
                    else
                        var_name = "z_bar" * string(d - nm)
                    end
                    
                    monomial_str *= "*" * var_name
                    if Avector[d] > 1
                        monomial_str *= "^" * string(Avector[d])
                    end
                end
            end

            for j in 1:nm
                coeff = Cp[p].fr[j, c]
                if abs(coeff) > 1e-20
                    coeff_str = "(" * string(round(real(coeff), digits=15)) * " + " * string(round(imag(coeff), digits=15)) * "im)"
                    rdyn_z[j] *= " + " * coeff_str * monomial_str
                end
            end
        end
    end

    ofile = open("./Output/equations_complex.txt", "w")
    for i in 1:nm
        write(ofile, rdyn_z[i] * ";\n")
    end
    close(ofile)

    println("Write document success ./Output/equations_complex.txt")
    println("Ending of writing equation!")
    println("========================================")
end

"""
write_rdyn_with_aliased_conjugates!(info::DPIMParam, Cp::Vector{Parametrisation})
"""
function write_rdyn_with_aliased_conjugates!(info::DPIMParam, Cp::Vector{Parametrisation})
    println("========================================")
    println("Begin write Reduced dynamic equations!")
    println("norm form = ", info.nrom )
    nm = info.nm
    nz = info.nz

    rdyn_z = ["" for i in 1:nz]
    for i in 1:nz
        rdyn_z[i] = "z"*string(i)*"' = "
    end

    for p in 1:info.max_order
        for c in 1:Cp[p].nc
            Avector = Cp[p].Avector[c]

            for j in 1:nm
                coeff = Cp[p].fr[j, c]

                if abs(coeff) < 1e-20
                    continue
                end

                monomial_str = ""
                for d in 1:info.nrom
                    if Avector[d] != 0
                        var_name = "z"*string(d)
                        monomial_str *= "*"*var_name
                        if Avector[d] > 1
                            monomial_str *= "^"*string(Avector[d])
                        end
                    end
                end
                coeff_str = "("*string(round(real(coeff), digits=15))*" + "*string(round(imag(coeff), digits=15))*"im)"
                rdyn_z[j] *= " + "*coeff_str*monomial_str

                conj_coeff = conj(coeff)
                conj_coeff_str = "("*string(round(real(conj_coeff), digits=15))*" + "*string(round(imag(conj_coeff), digits=15))*"im)"

                monomial_conj_str = ""
                for d in 1:info.nrom
                    if Avector[d] != 0
                        var_name = ""
                        if d <= nm
                            var_name = "z"*string(d + nm)
                        else
                            var_name = "z"*string(d - nm)
                        end
                        monomial_conj_str *= "*"*var_name
                        if Avector[d] > 1
                            monomial_conj_str *= "^"*string(Avector[d])
                        end
                    end
                end
                rdyn_z[j+nm] *= " + "*conj_coeff_str*monomial_conj_str
            end
        end
    end
    for i = 1 : info.nz
        println(rdyn_z[i])
    end
    filepath = "./Output/equations_aliased.txt"
    ofile = open(filepath, "w")
    for i in 1:nz
        write(ofile, rdyn_z[i]*";\n")
    end
    close(ofile)

    println("Write document success: ", filepath)
    println("Ending of writing equation!")
    println("========================================")
end

function SaveData(mappings,mappings_vel,mappings_modal,mappings_modal_vel,Avector,fdyn,howmany,ForceNode)
  matwrite("./Output/mappings.mat", Dict("mappings"=>mappings))
  matwrite("./Output/mappings_vel.mat", Dict("mappings_vel"=>mappings_vel))
  matwrite("./Output/mappings_modal.mat", Dict("mappings_modal"=>mappings_modal))
  matwrite("./Output/mappings_modal_vel.mat", Dict("mappings_modal_vel"=>mappings_modal_vel))
  matwrite("./Output/Avector.mat", Dict("Avector"=>Avector))
  matwrite("./Output/fdyn.mat", Dict("fdyn"=>fdyn))
  matwrite("./Output/howmany.mat", Dict("howmany"=>howmany))
  matwrite("./Output/ForceNode.mat", Dict("ForceNode"=>ForceNode))
  
end


  