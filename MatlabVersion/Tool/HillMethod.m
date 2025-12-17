function  [StableIndex, FloquetExponents] = HillMethod(system, dRdX)
% ref: The Harmonic Balance Method for Bifurcation Analysis of Large-Scale
% Nonlinear Mechanical Systems 
    sdof = system.sdof;
    B1 = -dRdX;
    B2 = kron(system.I, system.B);
    % =============================================
    B = B2 \ B1;
    % 
    [~, D] = eigs(B,size(B,1),'sm');
    Lambda = diag(D);
    % 
    image_Lambda = imag(Lambda);
    real_Lambda = real(Lambda);

    [~, Index] = sort(abs(image_Lambda));
    real_Lambda = real_Lambda(Index);
    FloquetExponents = real_Lambda(1:sdof);

    if all(FloquetExponents < 0)
        StableIndex = 1; 
    else
        StableIndex = 0; 
    end
end