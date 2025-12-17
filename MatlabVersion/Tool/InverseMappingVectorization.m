function [phys_disp, phys_vel] = InverseMappingVectorization(z, mappings_disp, mappings_vel, Avector)
% InverseMapping transforms reduced model coordinates to the real physical space
% via inverse mapping (vectorized implementation that handles negative z values).
%
% Inputs:
%   z              : Reduced model coordinate matrix, size [n_reduced x nTime]
%   mappings_disp  : Displacement inverse mapping coefficients, size [n_poly x n_nodes x n_dof]
%   mappings_vel   : Velocity inverse mapping coefficients, size [n_poly x n_nodes x n_dof]
%   Avector        : Polynomial power exponent matrix, size [n_poly x n_reduced]
%
% Outputs:
%   phys_disp      : Physical displacements, size [n_nodes x n_dof x nTime]
%   phys_vel       : Physical velocities, size [n_nodes x n_dof x nTime]

    % Get the dimensions of the variables
    [nPoly, nNodes, nDOF] = size(mappings_disp);
    [nPoly2, nReduced] = size(Avector);
    if nPoly ~= nPoly2
        error('Mismatch in the number of polynomial terms: The size of the first dimension of mappings_disp and Avector does not match.');
    end
    nTime = size(z, 2);
    
    % Calculate polynomial factors (handling negative values)
    abs_z = abs(z);       
    sign_z = sign(z); 
    
    polyFactors_abs = exp( Avector * log(abs_z) );  
    
    % Calculate the contribution from the signs
    polyFactors_sign = ones(nPoly, nTime);
    for j = 1:nReduced
        sign_part_j = repmat(sign_z(j,:), nPoly, 1) .^ repmat(Avector(:,j), 1, nTime);
        polyFactors_sign = polyFactors_sign .* sign_part_j;
    end
    
    % Final polynomial factors
    polyFactors = polyFactors_abs .* polyFactors_sign;  % size [nPoly x nTime]
    
    % Calculate the physical space data using matrix multiplication
    M_disp = reshape(mappings_disp, nPoly, nNodes*nDOF);
    M_vel  = reshape(mappings_vel,  nPoly, nNodes*nDOF);
    
    % Calculate physical displacements and velocities for all time steps using matrix multiplication.
    phys_disp_vec = M_disp' * polyFactors;
    phys_vel_vec  = M_vel'  * polyFactors;
    
    phys_disp = reshape(phys_disp_vec, nNodes, nDOF, nTime);
    phys_vel  = reshape(phys_vel_vec,  nNodes, nDOF, nTime);
end