function [Cp, TargetFreq] = ComputesParametrization(Mesh,SSMParams,MatParams)
    tic
    disp(["Initializing directories"])
    %
    folderName = 'Output';
    MakeFolder(folderName);
    disp(["============================================="])
    disp(["Assemblying M K"])
    [K,M] = AssembleKM(Mesh, MatParams);
    varInfo = whos('K');
    varMemoryMB = varInfo.bytes / 2^20;
    fprintf('Variable "K" occupies %.2f MB of memory.\n', varMemoryMB);
    varInfo = whos('M');
    varMemoryMB = varInfo.bytes / 2^20;
    fprintf('Variable "M" occupies %.2f MB of memory.\n', varMemoryMB);    
    disp(["============================================="])
    disp(["Computing eigenvalues"])
    [Phi, Lambda] = eigs(K(Mesh.Freedof,Mesh.Freedof),M(Mesh.Freedof,Mesh.Freedof),SSMParams.ComputeMode,'sm');
    Lambda = diag(Lambda);
    Phi = real(Phi);
    freq = sqrt(real(Lambda)); 
    TargetFreq = freq(SSMParams.MasterMode);
    disp([num2str(freq')])
    disp(["============================================="])
    C = MatParams.alpha*M+MatParams.beta*K; 
    % 
    PhiFull = zeros(Mesh.Sdof,SSMParams.ComputeMode);
    PhiFull(Mesh.Freedof,:) = Phi;
    for i = 1 : SSMParams.ComputeMode
        PhiSelectFull = PhiFull(:, i);
        PhiSelectFull = reshape(PhiSelectFull, 6, []);
        if sum(PhiSelectFull(1:3,:),"all")<0 || abs(min(PhiSelectFull(1:3,:),[],'all'))>abs(max(PhiSelectFull(1:3,:),[],'all'))
           Phi(:,i) = -Phi(:,i);
        end
    end
    disp(["Init Parametrisation"])
    Cp = InitParametrization(SSMParams);
    varInfo = whos('Cp');
    varMemoryMB = varInfo.bytes / 2^20;
    fprintf('Variable "Cp" occupies %.2f MB of memory.\n', varMemoryMB);  
    disp(["============================================="])
    % 
    [Cp,Sol_FULL,Rhs_FULL,Mat_FULL,Tri,AY,XTA] = SSMParamFirstOrder(SSMParams,MatParams,Mesh.ForceVector, Cp,freq,Phi,...
        K(Mesh.Freedof,Mesh.Freedof),M(Mesh.Freedof,Mesh.Freedof),C(Mesh.Freedof,Mesh.Freedof));
    % 
    [Cp] = SSMParamHigherOrder(Mesh,SSMParams,MatParams,Cp,freq,Phi,...
        K(Mesh.Freedof,Mesh.Freedof),M(Mesh.Freedof,Mesh.Freedof),C(Mesh.Freedof,Mesh.Freedof),...
        Sol_FULL,Rhs_FULL,Mat_FULL,Tri,AY,XTA);
    disp(["============================================="])
    toc
    disp(['Compute Paramters of invariant mainfold finish!'])
    disp(["============================================="])
end