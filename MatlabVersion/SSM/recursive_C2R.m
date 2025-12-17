function Cp = recursive_C2R(Ivec,p,pos,posinit,Avec,coeff,Cp,SSMParams)
    nzhalf=SSMParams.nz/2;
    nzfhalf=SSMParams.nForce;
    if pos==(p+1)
        pos1 = find(all(Cp.Avector == repmat(Avec, size(Cp.Avector, 1), 1), 2));
        Cp.Wr(:,pos1)= Cp.Wr(:,pos1) + coeff*Cp.W(:,posinit);
        Cp.fr(:,pos1)= Cp.fr(:,pos1) + coeff*Cp.f(:,posinit);
    else
        Avec1=Avec;
        Avec2=Avec;
        iz=Ivec(pos);
        if iz <= nzhalf
            coeff1=0.5*coeff;
            Avec1(iz) = Avec1(iz) + 1;
            coeff2 = -0.5*1i*coeff;
            Avec2(iz+nzhalf) = Avec2(iz+nzhalf) + 1;
        elseif iz <= SSMParams.nz
            coeff1=0.5*coeff;
            Avec1(iz-nzhalf) = Avec1(iz-nzhalf) + 1;
            coeff2 = 0.5*1i*coeff;
            Avec2(iz) = Avec2(iz) + 1;
        elseif iz <= SSMParams.nz+nzfhalf
            coeff1=1i*coeff;
            Avec1(iz) = Avec1(iz) + 1;
            coeff2=coeff;
            Avec2(iz+nzfhalf) = Avec2(iz+nzfhalf) + 1;
        else
            coeff1=-1i*coeff;
            Avec1(iz-nzfhalf) = Avec1(iz-nzfhalf) + 1;
            coeff2=coeff;
            Avec2(iz) = Avec2(iz) + 1;
        end
        pos = pos + 1;
        Cp = recursive_C2R(Ivec,p,pos,posinit,Avec1,coeff1,Cp,SSMParams);
        Cp = recursive_C2R(Ivec,p,pos,posinit,Avec2,coeff2,Cp,SSMParams);
    end
end