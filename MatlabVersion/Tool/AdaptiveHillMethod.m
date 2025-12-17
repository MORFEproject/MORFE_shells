function [StableIndex, FloquetExponents] = AdaptiveHillMethod(reducedSystem, dRdX, selectedHarmonics)

    sdof = reducedSystem.sdof;
    
    B1 = -dRdX;
    B2 = kron(reducedSystem.I, reducedSystem.B);
    B = B2 \ B1; 
    
    % 
    [~, D] = eigs(B, size(B,1), 'sm');
    Lambda = diag(D);
    
    %
    image_Lambda = imag(Lambda);
    real_Lambda = real(Lambda);
    
    % Moore filtre
    [~, Index] = sort(abs(image_Lambda));
    real_Lambda = real_Lambda(Index);
    
    %
    FloquetExponents = real_Lambda(1:sdof);

    % stable jugement
    if all(FloquetExponents < 0)
        StableIndex = 1; % stable
    else
        StableIndex = 0; % unstable
    end
end