function K = HBStiffnessMatrixLinear(Mext,Cext,Kext,omega,ScaleFactor1, ScaleFactor2)
    K = ScaleFactor1*ScaleFactor2^2*omega^2*Mext + ...
        ScaleFactor1*ScaleFactor2*omega*Cext + ScaleFactor1*Kext;
end