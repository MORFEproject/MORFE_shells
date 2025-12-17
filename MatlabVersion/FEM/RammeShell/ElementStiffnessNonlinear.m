function [kT,kNL] = ElementStiffnessNonlinear(CurrentElement, u, Mesh)
    kT = zeros(48, 48);
    kNL = zeros(48, 48);    
    [~, ~, nPtGaus] = GassPoint('2*2*2');

    for i = 1 : nPtGaus      
        G = Mesh.PreCalculateParam.G{CurrentElement, i};
        B0 = Mesh.PreCalculateParam.B0{CurrentElement, i};
        BNL = StrainMatrixNonlinear(G, u);
        det_a = Mesh.PreCalculateParam.det_a(CurrentElement, i);
        D = Mesh.PreCalculateParam.D{CurrentElement, i};
        sigma = D * (B0 + 0.5 * BNL) * u;
        S_extend = StressMatrix(sigma);
        kNL = kNL + 0.5*B0'*D*BNL*det_a + BNL'*D*B0*det_a + 0.5*BNL'*D*BNL*det_a;
        kT = kT + B0'*D*BNL*det_a + BNL'*D*B0*det_a + BNL'*D*BNL*det_a + G'*S_extend*G*det_a;
    end

end