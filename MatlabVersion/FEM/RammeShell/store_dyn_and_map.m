function [mappings,mappings_vel,mappings_modal,mappings_modal_vel,Avector,fdyn] = ...
    store_dyn_and_map(Mesh,MatParams,SSMParams,Cp,howmany)
    tic
    disp(["Begin postprocess"])
    mappings = zeros(howmany, Mesh.NodeNum, 6);
    mappings_vel = zeros(howmany, Mesh.NodeNum, 6);
    mappings_modal = zeros(howmany,SSMParams.ComputeMode);
    mappings_modal_vel = zeros(howmany,SSMParams.ComputeMode);
    Avector = zeros(howmany,Cp{1}.nc);
    fdyn=zeros(howmany,SSMParams.nm*2);
    disp(["Assemblying M K"])
    [K,M] = AssembleKM(Mesh, MatParams);
    [Phi, Lambda] = eigs(K(Mesh.Freedof,Mesh.Freedof),M(Mesh.Freedof,Mesh.Freedof),SSMParams.ComputeMode,'sm');
    Lambda = diag(Lambda);
    Phi = real(Phi);
    freq = sqrt(real(Lambda)); 
    disp(['Computing eigenvalues'])
    disp([num2str(freq')])
    PhiFull = zeros(Mesh.Sdof,SSMParams.ComputeMode);
    PhiFull(Mesh.Freedof,:) = Phi;
    for i = 1 : SSMParams.ComputeMode
        PhiSelectFull = PhiFull(:, i);
        PhiSelectFull = reshape(PhiSelectFull, 6, []);
        if sum(PhiSelectFull(1:3,:),"all")<0 || abs(min(PhiSelectFull(1:3,:),[],'all'))>abs(max(PhiSelectFull(1:3,:),[],'all'))
           Phi(:,i) = -Phi(:,i);
        end
    end 
    index=1;
    for p = 1 : SSMParams.max_order
        for i = 1:Cp{p}.nc
            Avector(index,:)=Cp{p}.Avector(i,:);
            index=index+1;
        end
    end
    index=1;
    for p = 1:SSMParams.max_order
        for i = 1:Cp{p}.nc
            for j = 1:SSMParams.nm
                fdyn(index,j)=2*real(Cp{p}.fr(j,i));
            end
            for j = 1:SSMParams.nm
                fdyn(index,j+SSMParams.nm)=-2*imag(Cp{p}.fr(j,i));
            end
            index=index+1;
        end
    end
    index=1;
    Dof = FindDof(Mesh.Sdof,Mesh.Freedof);
    for p = 1:SSMParams.max_order 
        for i = 1:Cp{p}.nc 
            %
            for inode = 1:Mesh.NodeNum
                for idof=1:6
                    dof = Dof(6*(inode-1)+idof);
                    if dof>0 
                        mappings(index,inode,idof)=real(Cp{p}.Wr(dof,i)); % displacemenet
                        mappings_vel(index,inode,idof)=real(Cp{p}.Wr(SSMParams.nK+dof,i)); % velocity
                    end
                end
            end
            %
            for imode=1:SSMParams.ComputeMode
                mappings_modal(index,imode)=Phi(:,imode)'*...
                    M(Mesh.Freedof,Mesh.Freedof)*...
                    real(Cp{p}.Wr(1:SSMParams.nK,i)); 
                mappings_modal_vel(index,imode)=Phi(:,imode)'*...
                    M(Mesh.Freedof,Mesh.Freedof)*...
                    real(Cp{p}.Wr(SSMParams.nK+1:2*SSMParams.nK,i)); 
            end
            index=index+1;
        end
    end
    toc
    disp(["End postprocess!"])
    disp(["============================================="])
end