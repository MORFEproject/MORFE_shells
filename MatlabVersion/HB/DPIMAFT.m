function [dfdx, f, U] = DPIMAFT(system, HBParam, FourierInverseTerms)
    Q   = system.Q;
    NT  = HBParam.N_T;
    sdof = system.sdof;
    
    % 
    Gamma   = FourierInverseTerms.Gamma;
    Gamma_p = FourierInverseTerms.Gamma_p;
    
    U = Gamma * Q;  
    U = reshape(U, [sdof, NT]);  
    
    f = zeros(sdof, NT);
    dfdxMatrix = zeros(sdof*NT, sdof*NT);
    
    for i = 1:NT
        [fnlt, KT] = AssembleNonlinear_mex(U(:,i), system.eqStructNumeric, i, NT);
        f(:,i) = fnlt;
        idx = (i-1)*sdof + (1:sdof);
        dfdxMatrix(idx, idx) = KT;
    end
    
    % 
    f = reshape(f, [], 1);
    f = Gamma_p * f;
    
    %
    temp = Gamma_p * dfdxMatrix;
    dfdx = temp * Gamma;
end
