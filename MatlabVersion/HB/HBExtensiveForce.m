function Fext = HBExtensiveForce(F, NH)
    sdof = length(F);
    Fext = sparse((2*NH+1)*sdof,1);
    for i = 1 : sdof
        Fext(sdof+i) = F(i);
    end
end