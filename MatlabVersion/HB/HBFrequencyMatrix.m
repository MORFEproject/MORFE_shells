function Rmc = HBFrequencyMatrix(Mext,Cext,omega,ScaleFactor1, ScaleFactor2)
    Rmc = 2*ScaleFactor1*ScaleFactor2^2*omega*Mext + ...
        ScaleFactor1*ScaleFactor2*Cext;
end