function [Q, Omega, Y_cont, loadStepConverged,StableSystem] = DPIMArclengthMethodImproved(system, HBParam)
tic
    Q = sparse((2*HBParam.N_H+1)*system.sdof,HBParam.Iter_max);
    Omega  = sparse(1,HBParam.Iter_max);
    Omega(1) = system.omega*system.ScaleFactor2;
    Q(:, 1) = system.Q*system.ScaleFactor1;
    Y_cont   = sparse(system.sdof,HBParam.Iter_max);
    psi = 1.0; % arc-length param if psi = 1, spherical arc length method

    % 
    direction = sign(HBParam.freq_end - HBParam.freq_start);

    % initial displacement
    dispPrev = zeros((2*HBParam.N_H+1)*system.sdof,1);
    dispPrev2 = zeros((2*HBParam.N_H+1)*system.sdof,1);

    % arc-length increment
    Ds = HBParam.loadincr;
    DsPrev = Ds;
    DsMax = Ds;
    DsMin = Ds;

    % load param increment 
    loadfactorPrev2 = system.omega*1;
    loadfactorPrev  = system.omega*1;
    
    % converge judgement
    converged = false;
    
    loadStepConverged = 0;
    for loadStep = 1 : HBParam.Iter_max
        % predict step
        if loadStep > 1
            DsFactor1 = Ds/DsPrev;
            system.Q = (1.0+DsFactor1)*dispPrev - DsFactor1*dispPrev2;
            system.omega = (1.0+DsFactor1)*loadfactorPrev - DsFactor1*loadfactorPrev2; 
        end 
        if loadStep > 1 && converged
            dRdXOld = dRdXNew;
            [StableOld, FloquetExponentsOld] = HillMethod(system, dRdXOld);
        end
        Du = system.Q - dispPrev; 
        Dl = system.omega - loadfactorPrev; 
        convergedPrev = converged;
        converged = false;   

        % correct steps
        for iter = 1 : HBParam.k_max
            % AFT
            [dfdx, f, U] = DPIMAFT(system, HBParam, system.FourierInverseTerms);
            Kmc0 = system.ScaleFactor1*system.ScaleFactor2*...
                system.omega*kron(system.delta1, system.B) +...
                system.ScaleFactor1*kron(system.I, system.A);
            K = Kmc0 - system.ScaleFactor1*dfdx;
            % 
            Rmc = system.ScaleFactor1*system.ScaleFactor2*...
                kron(system.delta1, system.B)* system.Q;
            R = -system.ScaleFactor1*system.ScaleFactor2*...
                system.omega*...
                kron(system.delta1, system.B)*system.Q +...
                system.ScaleFactor1*f;
            [converged, du, dl] = SolveArcLength(loadStep, iter, K, R, Rmc, Du, Dl, Ds, psi);
            if converged
                break;
            end    
            system.Q = system.Q + du;
            system.omega = system.omega + dl;
            Du = Du + du;
            Dl = Dl + dl;
        end

        if converged
            dRdXNew = K;
            [StableNew, FloquetExponentsNew] = HillMethod(system, dRdXNew);
            if(loadStep == 1)
                Ds = HBParam.loadincr;
                DsMax = Ds*10;
                DsMin = Ds/100000;
                FloquetExponentsOld = FloquetExponentsNew;
                dRdXOld = dRdXNew;
                StableOld = StableNew;
            end
            DsPrev = Ds;
            if convergedPrev
                Ds = min(max(2.0*Ds, DsMin), DsMax);
            end
            dispPrev2 = dispPrev;
            dispPrev = system.Q;
            loadfactorPrev2 = loadfactorPrev;
            loadfactorPrev = system.omega;
            loadStepConverged = loadStepConverged + 1;
            Q(:, loadStepConverged) = system.Q*system.ScaleFactor1;
            Omega(loadStepConverged) = system.omega*system.ScaleFactor2;
            Y_cont(:, loadStepConverged)  = max(U*system.ScaleFactor1,[],2); 

            % ======================== BifurcationType ========================
            BifurcationType = CheckBifurcationType(FloquetExponentsNew, FloquetExponentsOld, dRdXNew, dRdXOld);
            StableSystem{loadStepConverged}.StableIndex = StableOld;
            StableSystem{loadStepConverged}.BifurcationType = BifurcationType;
            % =============================================================
            disp(['============================='])
            disp(['Current step reach convergence!'])
            disp(['Current step:', num2str(loadStepConverged), ', Step size:', num2str(Ds)])
            disp(['Current Omega:', num2str(system.omega*system.ScaleFactor2)])
            
            % 
            if direction > 0
                if system.omega < 0 || system.omega*system.ScaleFactor2 > HBParam.freq_end
                    disp(['Reached bounded frequencies: ' num2str(system.omega*system.ScaleFactor2)])
                    break;
                end
            else % backward continuation: freq_end < freq_start
                if system.omega*system.ScaleFactor2 < HBParam.freq_end
                    disp(['Reached bounded frequencies: ' num2str(system.omega*system.ScaleFactor2)])
                    break;
                end
            end
        else
            if convergedPrev
                Ds = max(Ds*0.5, DsMin);
            else
                Ds = max(Ds*0.25, DsMin);
            end
        end
    end
toc
end
