function [converged, du, dl, du1] = SolveArcLength(timeStep, iter, K, R, Rmc, Du, Dl, ds, psi)
    if  timeStep > 1
        A = -(Du'*Du + psi*Dl*Dl - ds*ds);
        a = 2.0*Du';
        b = 2.0*psi*Dl;
    else
        A = 0.0;
        a = 0.0*Du';
        b = 1.0;
    end

    rNorm = norm(R,2); 
    rNorm = rNorm + norm(A); 

    fprintf(' rNorm : %5d ...  %12.6E \n', iter, rNorm);
    du = R*0.0;
    dl = 0.0;
    converged = false;

    if rNorm < 1.0e-8
       converged = true;
       return;
    end

    du1 = K\Rmc;
    du2 = K\R;
    dl = (a*du2 - A)/(a*du1-b);
    du = du2 - dl*du1;
end