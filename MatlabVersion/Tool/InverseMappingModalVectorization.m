function [modal_disp, modal_vel] = InverseMappingModalVectorization(z, mappings_modal_disp, mappings_modal_vel, Avector)
    
    [n_poly, n_modal] = size(mappings_modal_disp);
    [n_poly2, n_reduced] = size(Avector);
    if n_poly ~= n_poly2
        error('Mismatch in the number of polynomial terms: The size of the first dimension of mappings_disp and Avector does not match.');
    end
    
    nTime = size(z, 2);
    
    abs_z = abs(z);      
    sign_z = sign(z);     
   
    polyFactors_abs = exp( Avector * log(abs_z) );  
    polyFactors_sign = ones(n_poly, nTime);
    for j = 1:n_reduced
        polyFactors_sign = polyFactors_sign .* ( repmat(sign_z(j,:), n_poly, 1) .^ repmat(Avector(:,j), 1, nTime) );
    end
    
    polyFactors = polyFactors_abs .* polyFactors_sign;  
    
    modal_disp = mappings_modal_disp' * polyFactors;
    modal_vel  = mappings_modal_vel'  * polyFactors;
    
end
