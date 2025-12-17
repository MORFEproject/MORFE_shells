function [Cp,Sol,Rhs,Mat,Tri,AY,XTA] = SSMParamFirstOrderHalf(SSMParams,MatParams,ForceVector, Cp,freq,Phi,K,M,C)
    disp(["Order 1"])
    
    Rhs = zeros(SSMParams.nK+SSMParams.nz,1);
    Sol = zeros(SSMParams.nK+SSMParams.nz,1);
    Mat= zeros(SSMParams.nK+SSMParams.nz,SSMParams.nK+SSMParams.nz);

    omega0 = zeros(SSMParams. nm,1); % natural frequency
    elta0 = zeros(SSMParams.nm,1); % damping param
    AY = zeros(SSMParams.nA,SSMParams.nz);
    XTA = zeros(SSMParams.nz,SSMParams.nA);
    for i = 1 : SSMParams.nm
        omega0(i) = freq(SSMParams.MasterMode(i));
        elta0(i) = 0.5*(MatParams.alpha/omega0(i) + MatParams.beta*omega0(i));
        disp(['Damping param = ',num2str(elta0(i))])
        lambda1 = -elta0(i)*omega0(i) + omega0(i)*sqrt(complex(1-elta0(i)^2))*1i;
        lambda2 = -elta0(i)*omega0(i) - omega0(i)*sqrt(complex(1-elta0(i)^2))*1i;
        disp(['Damping param modified1 = ',num2str(lambda1)])
        disp(['Damping param modified2 = ',num2str(lambda2)])
        % 
        Cp{1}.f(i,i) = lambda1;% 
        Cp{1}.f(i+SSMParams.nm,i+SSMParams.nm) = lambda2;       
        Cp{1}.W(1:SSMParams.nK,i) = Phi(:,SSMParams.MasterMode(i));
        Cp{1}.W(SSMParams.nK+1:SSMParams.nA,i) = Phi(:,SSMParams.MasterMode(i)) * lambda1;
        Cp{1}.W(1:SSMParams.nK,i+SSMParams.nm) = Phi(:,SSMParams.MasterMode(i));
        Cp{1}.W(SSMParams.nK+1:SSMParams.nA,i+SSMParams.nm) = Phi(:,SSMParams.MasterMode(i)) * lambda2;
        
        AY(1:SSMParams.nK,i) = M*Cp{1}.W(1:SSMParams.nK,i);
        AY(1:SSMParams.nK,i+SSMParams.nm) = AY(1:SSMParams.nK,i);
        AY(SSMParams.nK+1:SSMParams.nA,i) = AY(1:SSMParams.nK,i)*lambda1;
        AY(SSMParams.nK+1:SSMParams.nA,i+SSMParams.nm) = AY(1:SSMParams.nK,i)*lambda2;
        XTA(i,1:SSMParams.nK) = -AY(1:SSMParams.nK,i)*lambda2/(lambda1-lambda2);
        XTA(i+SSMParams.nm,1:SSMParams.nK) = -AY(1:SSMParams.nK,i)*lambda1/(lambda2-lambda1);
        XTA(i,SSMParams.nK+1:SSMParams.nA) = AY(1:SSMParams.nK,i)/(lambda1-lambda2);
        XTA(i+SSMParams.nm,SSMParams.nK+1:SSMParams.nA) = AY(1:SSMParams.nK,i)/(lambda2-lambda1);
    end
    % external force
    Tri = zeros(SSMParams.nrom,1);
    for i = 1 : SSMParams.nz
        Tri(i)=Cp{1}.f(i,i); 
    end
    if SSMParams.nForce==1
        Tri(SSMParams.nz+1)=1i*omega0(SSMParams.Ffreq)*SSMParams.omega_mul;
        Cp{1}.f(SSMParams.nz+1,SSMParams.nz+1) = Tri(SSMParams.nz+1);
        Tri(SSMParams.nz+2)= - 1i*omega0(SSMParams.Ffreq)*SSMParams.omega_mul;
        Cp{1}.f(SSMParams.nz+2,SSMParams.nz+2) = Tri(SSMParams.nz+2);
    end
    for i = 1 : 2*SSMParams.nForce
       for j =  1 : size(SSMParams.Fmodes,1)
           Cp{1}.rhs(SSMParams.nK+1:SSMParams.nA,SSMParams.nz+i)= Cp{1}.rhs(SSMParams.nK+1:SSMParams.nA,SSMParams.nz+i)+...
           SSMParams.Fmult(j)*M*Phi(:,SSMParams.Fmodes(j));
       end
           [Sol,Rhs,Mat,Cp] = ...
           homological_HALF(Rhs,Mat,Tri,Cp,1,SSMParams.nz+i,M,K,C,AY,SSMParams,MatParams);
           % [Sol_FULL,Rhs_FULL,Mat_FULL,Cp{1}] = ...
           %     homological_Full(Rhs_FULL,Mat_FULL,Tri,Cp{1},SSMParams.nz+i,M,K,C,AY,XTA,SSMParams);
    end
end