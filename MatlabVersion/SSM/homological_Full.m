function [Sol,Rhs,Mat,Cp] = homological_Full(Rhs,Mat,Tri,Cp,Apos,M,K,C,AY,XTA,SSMParams)
    sigma=dot(Cp.Avector(Apos,:),Tri(1:SSMParams.nrom));
    resonant_modes = zeros(SSMParams.nrom,1);
    resonant_modes(1:end) = 0;
    if SSMParams.style == 'c'
        for i = 1 : SSMParams.nz
           lambda =  Tri(i);
           if (abs(imag(sigma-lambda))/abs(imag(lambda))<=1e-3)
               resonant_modes(i) = 1;
           end
        end
    elseif SSMParams.style == 'g'
        resonant_modes(:) = 1;
    elseif SSMParams.style == 'r'
        for i = 1 : SSMParams.nz
            lambda =  Tri(i);
            if (abs(imag(sigma-lambda))/abs(imag(lambda))<=1e-3)
                resonant_modes(i) = 1;
                if (i <= (SSMParams.nz/2))
                    resonant_modes(i+(SSMParams.nz/2)) = 1;
                else
                    resonant_modes(i-(SSMParams.nz/2)) = 1;
                end
            end
        end
    end
    Rhs = sparse(SSMParams.nMat,1);
    disp(['Current MonomialExponents = ',num2str(Cp.Avector(Apos,:)),';Resonant index = ',num2str(resonant_modes')])
    Rhs(1:SSMParams.nA) = Cp.rhs(:,Apos);
    Rhs(1:SSMParams.nK) = Rhs(1:SSMParams.nK) - M*Cp.Wf(1:SSMParams.nK,Apos);
    Rhs(SSMParams.nK+1:SSMParams.nA) = Rhs(SSMParams.nK+1:SSMParams.nA) - M*Cp.Wf(SSMParams.nK+1:SSMParams.nA,Apos);
    Mat = sparse(SSMParams.nMat,SSMParams.nMat);
    Mat(1:SSMParams.nK,1:SSMParams.nK) = sigma*M;
    Mat(SSMParams.nK+1:SSMParams.nA,SSMParams.nK+1:SSMParams.nA) = sigma*M+C;
    Mat(1:SSMParams.nK,SSMParams.nK+1:SSMParams.nA) = -M;
    Mat(SSMParams.nK+1:SSMParams.nA,1:SSMParams.nK) = K;
    
    for j = 1 : SSMParams.nz
       if  resonant_modes(j) 
           Mat(1:SSMParams.nA,SSMParams.nA+j)=AY(:,j);
           Mat(SSMParams.nA+j,1:SSMParams.nA)=XTA(j,:);
       else
           Mat(SSMParams.nA+j,SSMParams.nA+j)=1;
       end
    end

    if SSMParams.nMat > 100000
        setup.type = 'ilutp';
        setup.milu = 'row';
        setup.droptol = 1e-8;
        [L, U] = ilu(Mat, setup);
        Sol = qmr(Mat, Rhs, 1e-6, 8000, L, U);
    else
        Sol=Mat\Rhs; 
    end

    Cp.W(:,Apos)=Sol(1:SSMParams.nA);
    Cp.f(1:SSMParams.nz,Apos)=Sol(SSMParams.nA+1:SSMParams.nA+SSMParams.nz);
end