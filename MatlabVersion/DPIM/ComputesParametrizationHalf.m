function [Cp, TargetFreq] = ComputesParametrizationHalf(Mesh,DPIMParams,MatParams)
    tic
    disp(["Initializing directories"])
    % creat folder
    folderName = 'Output';
    MakeFolder(folderName);
    disp(["============================================="])
    disp(["Assemblying M K"])
    [K,M] = AssembleKM(Mesh, MatParams);
    % output memory of K, M
    varInfo = whos('K');
    varMemoryMB = varInfo.bytes / 2^20;
    fprintf('Variable "K" occupies %.2f MB of memory.\n', varMemoryMB);
    varInfo = whos('M');
    varMemoryMB = varInfo.bytes / 2^20;
    fprintf('Variable "M" occupies %.2f MB of memory.\n', varMemoryMB);    
    fprintf('Total degree of freedom is %d.\n', size(K,1)); 
    fprintf('Free degree of freedom is %d.\n', Mesh.FreedofNum); 
    disp(["============================================="])
    disp(["Computing eigenvalues"])
    [Phi, Lambda] = eigs(K(Mesh.Freedof,Mesh.Freedof),M(Mesh.Freedof,Mesh.Freedof),DPIMParams.ComputeMode,'sm');
    Lambda = diag(Lambda);
    Phi = real(Phi);
    freq = sqrt(real(Lambda)); 
    TargetFreq = freq(DPIMParams.MasterMode);
    disp([num2str(freq')])
    disp(["============================================="])
    C = MatParams.alpha*M+MatParams.beta*K; % damping matrix 
    % 
    PhiFull = zeros(Mesh.Sdof,DPIMParams.ComputeMode);
    PhiFull(Mesh.Freedof,:) = Phi;
    for i = 1 : DPIMParams.ComputeMode
        PhiSelectFull = PhiFull(:, i);
        PhiSelectFull = reshape(PhiSelectFull, 6, []);
        if sum(PhiSelectFull(1:3,:),"all")<0 || abs(min(PhiSelectFull(1:3,:),[],'all'))>abs(max(PhiSelectFull(1:3,:),[],'all'))
           Phi(:,i) = -Phi(:,i);
        end
    end
    % 
    disp(["Init Parametrisation"])
    Cp = InitParametrization(DPIMParams);
    varInfo = whos('Cp');
    varMemoryMB = varInfo.bytes / 2^20;
    fprintf('Variable "Cp" occupies %.2f MB of memory.\n', varMemoryMB);  
    disp(["============================================="])
    % first order of dpim
    [Cp,Sol_FULL,Rhs_FULL,Mat_FULL,Tri,AY,XTA] = DPIMParamFirstOrderHalf(DPIMParams,MatParams,Mesh.ForceVector, Cp,freq,Phi,...
        K(Mesh.Freedof,Mesh.Freedof),M(Mesh.Freedof,Mesh.Freedof),C(Mesh.Freedof,Mesh.Freedof));
    % p-th order of dpim
    [Cp] = DPIMParamHigherOrderHalf(Mesh,DPIMParams,MatParams,Cp,freq,Phi,...
        K(Mesh.Freedof,Mesh.Freedof),M(Mesh.Freedof,Mesh.Freedof),C(Mesh.Freedof,Mesh.Freedof),...
        Sol_FULL,Rhs_FULL,Mat_FULL,Tri,AY,XTA);
    disp(["============================================="])
    toc
    disp(['Compute Paramters of invariant mainfold finish!'])
    disp(["============================================="])
end