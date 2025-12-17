function [system, HBParam, Om_HB,Q_HB,Ualpmlitude,Falpmlitude,StableSystem] = ...
    RunHBMethod(NH, Oms, Ome, Arclength, ...
    ScaleButton, Platform, varargin)

    folderName = 'Output';
    currentFolder = pwd;
    targetFolderFullPath = fullfile(currentFolder, folderName);
    pathString = genpath(targetFolderFullPath);
    addpath(pathString);
    load Output\classified_eq_numeric.mat 
    InitModule('DPIM','HB');
    switch Platform
        case 'NLvib'
            %% ========== Solve and continue (the part of NLvib) ===================
            HBParam.N_H = NH;  % HB order
            HBParam.N_T = 4*HBParam.N_H+1; % time sampling
            HBParam.freq_start = Oms;
            HBParam.freq_end   = Ome; 
            % FFT param
            sdof = length(eqStructNumeric);  
            [delta1, I] = DPIMHBExtensiveMatrix(HBParam.N_H); % inverse FFT param
            FourierInverseTerms = DPIMFourierInverse(sdof, HBParam);
            % 
            B = eye(sdof, sdof);
            A = zeros(sdof, sdof);
            % ======================== Matrix A(linear stiffness) =========================
            for i = 1 : sdof
                tempStruct = eqStructNumeric(i);  
                linTerms = tempStruct.linear;
                for r = 1 : size(linTerms, 1)
                    coeff = linTerms(r, 1);
                    exponents = linTerms(r, 2:end); 
                    physExponents = exponents(1:sdof);
                    idx = find(physExponents == 1);
                    if length(idx) == 1 && all(physExponents(setdiff(1:sdof, idx)) == 0)
                        % x' = B^{-1}(f_ext - A x)
                        A(i, idx) = -coeff;
                    end
                end
            end
            % ==================== the first order of external force ================
            U = zeros(sdof, HBParam.N_T); 
            size1 = size(U,1);
            size2 = size(U,2);
            fext = zeros(size1, size2);
            for i = 1:HBParam.N_T
                [fnext, ~] = AssembleNonlinear_mex(U(:,i), eqStructNumeric, i, HBParam.N_T);
                fext(:,i) = fnext;
            end
            Falpmlitude = max(fext,[],2);
            fext = reshape(fext, [], 1);
            Fext = FourierInverseTerms.Gamma_p * fext;
            system.Fext = Fext;
            % ==============================================================
            KextLinear = Oms*kron(delta1, B) + kron(I, A);
            Q0 = KextLinear\Fext; % inital guess
            % Building system data
            system.sdof = sdof;
            system.I = I;
            system.delta1 = delta1;
            system.FourierInverseTerms = FourierInverseTerms;
            
            system.B = B;
            system.A = A;
            system.eqStructNumeric = eqStructNumeric;
            ds = Arclength;
            analysis_type = 'frf';
        
            if ScaleButton == 0
                Sopt = struct('Dscale',[1e-0*ones(size(Q0));1]);
            else
                Sopt = struct('Dscale',[1e-0*ones(size(Q0));Oms]);
            end
            [X_HB,Solinfo_HB,Sol,stable,FloquetExponents,Bifurcation] = solve_and_continue(Q0,...
                @(X) HB_residual(X,system,HBParam.N_H, HBParam.N_T,analysis_type),...
                Oms,Ome,ds,system,Sopt);
            Om_HB = X_HB(end,:);
            Q_HB =  X_HB(1:end-1,:);
            Ualpmlitude = zeros(sdof,size(Q_HB,2));
            for i = 1 : size(Q_HB,2)
                U = FourierInverseTerms.Gamma * Q_HB(:,i);
                U = reshape(U, [sdof, HBParam.N_T]);
                Ualpmlitude(:, i) = max(U, [], 2);
            end

            StableSystem = cell(length(stable));
            for i = 1 : length(stable)
                StableSystem{i}.StableIndex = stable(i);
                StableSystem{i}.BifurcationType = Bifurcation{i};
            end
        case 'AdaptiveArclength'
            %% ========== Solve and continue (ADHB) ===================
            HBParam.N_H = NH;  % HB order
            HBParam.N_T = 4*HBParam.N_H+1; % time sampling
            HBParam.freq_start = Oms;
            HBParam.freq_end   = Ome; 
            % FFT param
            sdof = length(eqStructNumeric);  
            [delta1, I] = DPIMHBExtensiveMatrix(HBParam.N_H); % inverse FFT param
            FourierInverseTerms = DPIMFourierInverse(sdof, HBParam);
            % 
            B = eye(sdof, sdof);
            A = zeros(sdof, sdof);
            % ======================== Matrix A(linear stiffness) =========================
            for i = 1 : sdof
                tempStruct = eqStructNumeric(i);  
                linTerms = tempStruct.linear;
                for r = 1 : size(linTerms, 1)
                    coeff = linTerms(r, 1);
                    exponents = linTerms(r, 2:end); 
                    physExponents = exponents(1:sdof);
                    idx = find(physExponents == 1);
                    if length(idx) == 1 && all(physExponents(setdiff(1:sdof, idx)) == 0)
                        % x' = B^{-1}(f_ext - A x)
                        A(i, idx) = -coeff;
                    end
                end
            end
            % ==================== the first order of external force ================
            U = zeros(sdof, HBParam.N_T); 
            size1 = size(U,1);
            size2 = size(U,2);
            fext = zeros(size1, size2);
            for i = 1:HBParam.N_T
                [fnext] = AssembleNonlinearForce(U(:,i), eqStructNumeric, i, HBParam.N_T);
                fext(:,i) = fnext;
            end
            Falpmlitude = max(fext,[],2);
            fext = reshape(fext, [], 1);
            Fext = FourierInverseTerms.Gamma_p * fext;
            system.Fext = Fext;
            % ==============================================================
            KextLinear = Oms*kron(delta1, B) + kron(I, A);
            Q0 = KextLinear\Fext; % inital guess
            % Building system data
            system.sdof = sdof;
            system.I = I;
            system.delta1 = delta1;
            system.FourierInverseTerms = FourierInverseTerms;
            
            system.B = B;
            system.A = A;
            system.eqStructNumeric = eqStructNumeric;
            HBParam.loadincr = Arclength;
            HBParam.Iter_max = 4000;
            HBParam.k_max = 15;
            if ScaleButton == 0
                ScaleFactor1 = 1; 
                ScaleFactor2 = 1;
            else
                ScaleFactor1 = 1; 
                ScaleFactor2 = Ome;
            end
            system.ScaleFactor1 = ScaleFactor1;
            system.ScaleFactor2 = ScaleFactor2;
            system.Q = Q0/system.ScaleFactor1;
            system.omega = Oms/system.ScaleFactor2;
            [Q_HB, Om_HB, Ualpmlitude, IterationNum, StableSystem] = DPIMArclengthMethodImproved(system, HBParam);
            % [Q_HB, Om_HB, Ualpmlitude, IterationNum, StableSystem] = DPIMArclengthMethod(system, HBParam);
            Om_HB = Om_HB(1:IterationNum);
            Ualpmlitude = Ualpmlitude(:, 1:IterationNum);
            Q_HB = Q_HB(:, 1:IterationNum);
    end
end