function Fe = Integrate_H(PreCalculateParam,u1,u2,u3,MatParams,CurrentElementNum)
    Fe = zeros(48,1);
    [~, ~, PointNum] = GassPoint('2*2*2'); 
    kau_nl = zeros(4, 48);
    kua_nl = zeros(48, 4);
    for i = 1 :  PointNum
        BL1 = StrainMatrixNonlinear(PreCalculateParam.G{CurrentElementNum, i}, u1);
        BL2 = StrainMatrixNonlinear(PreCalculateParam.G{CurrentElementNum, i}, u2);      
        BL3 = StrainMatrixNonlinear(PreCalculateParam.G{CurrentElementNum, i}, u3);
        D = PreCalculateParam.D{CurrentElementNum, i};
        sigma12 = 0.5*D*BL1*u2;sigma21 = 0.5*D*BL2*u1;
        sigma13 = 0.5*D*BL1*u3;sigma31 = 0.5*D*BL3*u1;
        sigma23 = 0.5*D*BL2*u3;sigma32 = 0.5*D*BL3*u2; 
        kau_nl = kau_nl + PreCalculateParam.Ba{CurrentElementNum, i}'*...
            PreCalculateParam.D{CurrentElementNum, i}*...
            BL2*...
            PreCalculateParam.det_a(CurrentElementNum, i);
        kua_nl = kua_nl + BL1'*...
            PreCalculateParam.D{CurrentElementNum, i}*...
            PreCalculateParam.Ba{CurrentElementNum, i}*...
            PreCalculateParam.det_a(CurrentElementNum, i);
        Fe = Fe + (1/6)*(BL1'*sigma23+BL1'*sigma32+...
                         BL2'*sigma13+BL2'*sigma31+...
                         BL3'*sigma12+BL3'*sigma21)*PreCalculateParam.det_a(CurrentElementNum, i);
    end
    Fe = Fe - 0.5*kua_nl*PreCalculateParam.invKaa{CurrentElementNum, 1}*kau_nl*u3;
end