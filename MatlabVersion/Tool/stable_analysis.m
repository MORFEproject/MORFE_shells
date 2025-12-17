function [StableIndex,FloquetExponent,BifurcationType]= stable_analysis(DRDXold,DRDXnew,system,stable, FloquetExponents,istep)
    B = system.B;
    A = system.A;
    
    delta1 = system.delta1;
    I = system.I;
    sdof=length(B);
    dfdxnew=DRDXnew(:,1:end-1);
    dfdxold=DRDXold(:,1:end-1);
    B1new=-dfdxnew;
    B1old=-dfdxold;
    B2=kron(I, B);
    Binew=B2\B1new;
    Biold=B2\B1old;
    [~, Dnew] = eigs(Binew,length(Binew),'sm');
    
    Lambdanew = diag(Dnew);
    % 
    image_Lambdanew = imag(Lambdanew);
    real_Lambdanew = real(Lambdanew);
    % 
    [image_Lambdanew, Indexnew] = sort(abs(image_Lambdanew));
    real_Lambdanew = real_Lambdanew(Indexnew);
    if istep==2
        [~, Dold] = eigs(Biold,length(Biold),'sm');
        Lambdaold = diag(Dold);
        image_Lambdaold = imag(Lambdaold);
        real_Lambdaold = real(Lambdaold);
        [image_Lambdaold, Indexold] = sort(abs(image_Lambdaold));
        FloquetExponents=[FloquetExponents,real_Lambdaold(1:sdof)];
    end
    
    FloquetExponent = [FloquetExponents,real_Lambdanew(1:sdof)];
    
    if all(real_Lambdanew(1:sdof) < 0)
        if istep==2
            stable=[stable,1];
        end
        StableIndex = [stable;1]; 
    else
        if istep==2
            stable=[stable,0];
        end
        StableIndex = [stable;0];
    end
    
    %% bifurcation
    phi_SN_New = DetectionSN(dfdxnew);
    phi_SN_Old = DetectionSN(dfdxold);
    signSN = sign(phi_SN_New * phi_SN_Old);
    
    phi_NS_New = DetectionNS(FloquetExponent(:,istep));
    phi_NS_Old = DetectionNS(FloquetExponent(:,istep-1));
    signNS = sign(phi_NS_New * phi_NS_Old);
    
    if signSN <= 0
        BifurcationType = 'Saddle-node bifurcation';
    elseif signNS <=0
        BifurcationType = 'Neimark-Sacker bifurcation';
    else 
        BifurcationType = 'No Bifurcation';
    end
    
    % SN
    function phi_SN = DetectionSN(dRdX)
        phi_SN = det(dRdX);
    end
    
    % NS
    function phi_NS = DetectionNS(Floquet)
        n2 = length(Floquet); % n2 = 2n
        idxPairs = nchoosek(1:n2, 2);
        sums = Floquet(idxPairs(:,1)) + Floquet(idxPairs(:,2));
        phi_NS = prod(sums);
    end
    zxc=1;
end

