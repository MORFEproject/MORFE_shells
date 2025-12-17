function [A, Omega, Y_cont, loadStepConverged] = ArclengthMethod(system, HBParam)
%%
% Chennakesava Kadapa,
% A simple extrapolated predictor for overcoming the starting and tracking issues in the arc-length method for nonlinear structural mechanics,
% Engineering Structures,
% Volume 234,
% 2021,
% 111755,
% ISSN 0141-0296,
% https://doi.org/10.1016/j.engstruct.2020.111755.
%% ================ main code ========================
tic
    A = sparse((2*HBParam.N_H+1)*length(system.Mesh.Freedof),HBParam.Iter_max);
    Omega  = sparse(1,HBParam.Iter_max);
    Omega(1) = system.omega*system.ScaleFactor2;
    A(:, 1) = system.A*system.ScaleFactor1;
    Y_cont   = sparse(length(system.Mesh.Freedof),HBParam.Iter_max);
    FourierInverseTerms = FourierInverse(system, HBParam);
    psi = 1.0; % arc-length param if psi = 1, spherical arc length method

    % initial displacement
    dispPrev = zeros((2*HBParam.N_H+1)*length(system.Mesh.Freedof),1);
    dispPrev2 = zeros((2*HBParam.N_H+1)*length(system.Mesh.Freedof),1);

    % arc-lengence increment
    Ds = HBParam.loadincr;
    DsPrev = Ds;
    DsMax = Ds;
    DsMin = Ds;

    % load param increment (in this code, it means frequency increment)
    loadfactorPrev2 = system.omega*1;
    loadfactorPrev  = system.omega*1;
    
    % converge judgement
    converged = false;
    
    % solve direction
    isForwardSolve = HBParam.freq_start < HBParam.freq_end;
    
    loadStepConverged = 0;
    for loadStep = 1 : HBParam.Iter_max
        % predict step
        if  loadStep > 1
            DsFactor1 = Ds/DsPrev; 
            system.A = (1.0+DsFactor1)*dispPrev - DsFactor1*dispPrev2; 
            system.omega = (1.0+DsFactor1)*loadfactorPrev - DsFactor1*loadfactorPrev2; 
        end     
        Du = system.A - dispPrev;
        Dl = system.omega - loadfactorPrev; 
        convergedPrev = converged;
        converged = false;   

        % correct steps
        for iter = 1 : HBParam.k_max
            % AFT
            [dfdx, f, U] = AFTParallel(system, HBParam, FourierInverseTerms);
            % [dfdx, f, U] = AFT(system, HBParam, FourierInverseTerms);
            Kmc0 = HBStiffnessMatrixLinear(system.Mext,system.Cext,system.Kext, ...
                system.omega, system.ScaleFactor1, system.ScaleFactor2);
            K = Kmc0 + system.ScaleFactor1*dfdx;
            % 
            Rmc = HBFrequencyMatrix(system.Mext,system.Cext,system.omega, ...
                system.ScaleFactor1, system.ScaleFactor2);
            Rmc = Rmc * system.A;
            R = system.Fext - (Kmc0*system.A + system.ScaleFactor1*f);
            [converged, du, dl] = ...
                SolveArcLength(loadStep, iter, K, R, Rmc, Du, Dl, Ds, psi);
            if converged
                break;
            end    
            system.A = system.A + du;
            system.omega = system.omega + dl;
            Du = Du + du;
            Dl = Dl + dl;
        end

        if  converged
            if(loadStep == 1)
                % Ds = sqrt(Du'*Du + psi*Dl*Dl);
                Ds = HBParam.loadincr;
                DsMax = Ds*10;
                DsMin = Ds/10000;
            end
            DsPrev = Ds;
            if  convergedPrev
                Ds = min(max(2.0*Ds, DsMin), DsMax);
            end
            dispPrev2 = dispPrev;
            dispPrev = system.A;
            loadfactorPrev2 = loadfactorPrev;
            loadfactorPrev = system.omega;
            loadStepConverged = loadStepConverged + 1;
            A(:, loadStepConverged) = system.A*system.ScaleFactor1;
            Omega(loadStepConverged) = system.omega*system.ScaleFactor2;
            Y_cont(:, loadStepConverged)  = max(U*system.ScaleFactor1,[],2); 
            disp(['============================='])
            disp(['Current step reach converage!'])
            disp(['Current step:', num2str(loadStepConverged)])
            disp(['Current Omeage:', num2str(system.omega*system.ScaleFactor2)])
            
            % 
            if system.omega < 0 || (isForwardSolve && system.omega*system.ScaleFactor2 > HBParam.freq_end) || ...
               (~isForwardSolve && system.omega*system.ScaleFactor2 < HBParam.freq_end)
                disp(['reach bounded frequencies ' num2str(system.omega*system.ScaleFactor2)])
                break
            end       
        else
            if  convergedPrev
                Ds = max(Ds*0.5, DsMin);
            else
                Ds = max(Ds*0.25, DsMin);
            end
        end
    end
    A = A(:, 1:loadStepConverged);
toc
end