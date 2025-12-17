function Fe = Integrate_G(PreCalculateParam, u1, u2, MatParams, CurrentElementNum)
    Fe = zeros(48,1);
    [~, ~, PointNum] = GassPoint('2*2*2'); 
    kau_nl = zeros(4, 48);
    kua_nl = zeros(48, 4);
    for i = 1 :  PointNum
        BL1 = StrainMatrixNonlinear(PreCalculateParam.G{CurrentElementNum, i}, u1);
        BL2 = StrainMatrixNonlinear(PreCalculateParam.G{CurrentElementNum, i}, u2);
        B0 = PreCalculateParam.B0{CurrentElementNum, i};
        sigma1 = PreCalculateParam.D{CurrentElementNum, i} * B0 * u1;
        sigma2 = PreCalculateParam.D{CurrentElementNum, i} * B0 * u2;
        sigma12 = 0.5*PreCalculateParam.D{CurrentElementNum, i}*BL1*u2;
        sigma21 = 0.5*PreCalculateParam.D{CurrentElementNum, i}*BL2*u1;
        kau_nl = kau_nl + PreCalculateParam.Ba{CurrentElementNum, i}'*...
            PreCalculateParam.D{CurrentElementNum, i}*...
            BL1*...
            PreCalculateParam.det_a(CurrentElementNum, i); 
        kua_nl = kua_nl + BL1'*...
            PreCalculateParam.D{CurrentElementNum, i}*...
            PreCalculateParam.Ba{CurrentElementNum, i}*...
            PreCalculateParam.det_a(CurrentElementNum, i); 
        Fe = Fe + 0.5*(B0'*sigma12+B0'*sigma21+BL1'*sigma2+BL2'*sigma1)*...
            PreCalculateParam.det_a(CurrentElementNum, i);
    end
    Fe = Fe - 0.5*PreCalculateParam.K0a{CurrentElementNum, 1}*PreCalculateParam.invKaa{CurrentElementNum, 1}*kau_nl*u2 -...
        kua_nl*PreCalculateParam.invKaa{CurrentElementNum, 1}*PreCalculateParam.Ka0{CurrentElementNum, 1}*u2; 
end