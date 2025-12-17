function [Sol, Rhs, Mat, Cp] = homological_HALF(Rhs, Mat, Tri, Cp, p, Apos, M, K, C, AY, info, MatParams)

sigma = dot(Cp{p}.Avector(Apos, :), Tri(1:info.nrom));
resonant_modes = false(info.nrom, 1);
if strcmp(info.style, 'c')
    for i = 1:info.nz
        lambda = Tri(i);
        if (abs(imag(sigma - 1*lambda)) / abs(imag(1*lambda)) <= 1e-2)
            resonant_modes(i) = true;
        end
    end
elseif strcmp(info.style, 'g')
    resonant_modes(:) = true;
elseif strcmp(info.style, 'r')
    for i = 1:info.nz
        lambda = Tri(i);
        if abs(imag(sigma - lambda)) / abs(imag(lambda)) <= 1e-3
            resonant_modes(i) = true;
            if i <= (info.nz / 2)
                resonant_modes(i + info.nz / 2) = true;
            else
                resonant_modes(i - info.nz / 2) = true;
            end
        end
    end
end

disp(['Current MonomialExponents = ', num2str(Cp{p}.Avector(Apos,:)), '; Resonant index = ', num2str(resonant_modes')]);

Rhs = sparse(info.nMat, 1);
Rhs(1:info.nK) = Cp{p}.rhs(info.nK+1:info.nA, Apos) - M * Cp{p}.Wf(info.nK+1:info.nA, Apos) - (sigma * M + C) * Cp{p}.Wf(1:info.nK, Apos);

Mat = sparse(info.nMat, info.nMat);
Mat(1:info.nK, 1:info.nK) = (sigma^2 + MatParams.alpha * sigma) * M + (MatParams.beta * sigma + 1) * K;

for d = 1:info.nz
    Mat(info.nK+d, info.nK+d) = 1.0;
    if resonant_modes(d)
        lambda = Tri(d);
        Mat(info.nK+d, 1:info.nK) = (sigma - conj(lambda)) * AY(1:info.nK, d);
        Mat(1:info.nK, info.nK+d) = (sigma - conj(lambda)) * AY(1:info.nK, d);
        if d <= info.nm && resonant_modes(d + info.nm)
            Mat(info.nK+d, info.nK+d+info.nm) = 1.0;
        end
        if d > info.nm && resonant_modes(d - info.nm)
            Mat(info.nK+d, info.nK+d-info.nm) = 1.0;
        end
        Rhs(info.nK+d) = -AY(1:info.nK, d).' * Cp{p}.Wf(1:info.nK, Apos);
    end
end

Sol = Mat\Rhs;
Cp{p}.W(1:info.nK,Apos)=Sol(1:info.nK);
Cp{p}.W(info.nK+1:info.nA,Apos)=sigma*Sol(1:info.nK)+Cp{p}.Wf(1:info.nK,Apos);
Cp{p}.f(1:info.nz,Apos)=Sol(info.nK+1:info.nK+info.nz);
for d = 1:info.nz
    if resonant_modes(d) 
        Cp{p}.W(info.nK+1:info.nA,Apos) = Cp{p}.W(info.nK+1:info.nA,Apos) + Cp{p}.f(d,Apos)*Cp{1}.W(1:info.nK,d);
    end  
end

end
