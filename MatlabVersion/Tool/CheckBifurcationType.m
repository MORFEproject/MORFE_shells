function BifurcationType = CheckBifurcationType(FloquetExponentsNew, FloquetExponentsOld, dRdXNew, dRdXOld)
    phi_SN_New = DetectionSN(dRdXNew);
    phi_SN_Old = DetectionSN(dRdXOld);
    signSN = sign(phi_SN_New * phi_SN_Old);

    phi_NS_New = DetectionNS(FloquetExponentsNew);
    phi_NS_Old = DetectionNS(FloquetExponentsOld);
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
    function phi_NS = DetectionNS(FloquetExponents)
        n2 = length(FloquetExponents); % n2 = 2n
        idxPairs = nchoosek(1:n2, 2);
        sums = FloquetExponents(idxPairs(:,1)) + FloquetExponents(idxPairs(:,2));
        phi_NS = prod(sums);
    end

end
