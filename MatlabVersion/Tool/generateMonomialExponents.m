% This function is intended to replace the multiexponents method in the Julia Combinatorics library. 
% This method returns the exponents of each variable for a given polynomial. 
function exponents = generateMonomialExponents(nVars, degree)
    if nVars == 1
        exponents = [degree]; 
    else
        exponents = [];
        for d = degree:-1:0
            subExponents = generateMonomialExponents(nVars - 1, degree - d);
            exponents = [exponents; [repmat(d, size(subExponents, 1), 1), subExponents]];
        end
    end
end