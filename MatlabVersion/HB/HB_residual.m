function [R, dR] = HB_residual(X, system, H, NT, analysis_type,varargin)
    B = system.B;
    A = system.A;

    delta1 = system.delta1;
    I = system.I;
    
    sdof = system.sdof;
    % inverse discrete Fourier transform 
    FourierInverseTerms = system.FourierInverseTerms;
    eqStructNumeric = system.eqStructNumeric;



    % Handle analysis type
    if nargin<=3 || isempty(analysis_type)
        % Default analysis: frequency response
        analysis_type = 'frf';
    end
    switch lower(analysis_type)
        case {'frf','frequency response'}
            % FFT param
            Q = X(1:end-1);
            dQ = eye(length(X));
            dX = eye(length(X));
            Om = X(end);
            dOm = dX(end,:)	;
            dD_dalpha = 0*B; dalpha = zeros(1,length(X));
            % Scaling of dynamic force equilibrium
            if length(varargin)<2 || isempty(varargin{2})
                fscl = 1;
            else
                fscl = varargin{2};
            end
        otherwise
            error(['Unknown analysis type ' analysis.type '.']);
    end


    switch lower(analysis_type)
        case {'frf','frequency response'}
            % ================== AFT begin sin cos representation ======================
            U = FourierInverseTerms.Gamma * Q;
            U = reshape(U, [sdof, NT]);
            size1 = size(U,1);
            size2 = size(U,2);
            fnl = zeros(size1, size2);
            dfdxMatrix = zeros(sdof*NT, sdof*NT);
            for i = 1 : NT
                [fnlt, KT] = AssembleNonlinear_mex(U(:,i), eqStructNumeric, i, NT); 
                fnl(:,i) = fnlt;
                index_row_start = sdof*(i-1)+1;
                index_row_end = sdof*i;
                dfdxMatrix(index_row_start:index_row_end, index_row_start:index_row_end) = KT; 
            end
            fnl = reshape(fnl,[],1);
            fnl = FourierInverseTerms.Gamma_p * fnl;
            dfdx = FourierInverseTerms.Gamma_p * dfdxMatrix * FourierInverseTerms.Gamma;
            % ================== AFT end ======================
            R = Om*kron(delta1, B)*Q - fnl;
            dR = [Om*kron(delta1, B) + kron(I, A) - dfdx, kron(delta1, B)* Q];
        case {'nma','nonlinear modal analysis'}
            % ================== AFT begin sin cos representation ======================
            U = FourierInverseTerms.Gamma * Q;
            U = reshape(U, [sdof, NT]);
            size1 = size(U,1);
            size2 = size(U,2);
            fnl = zeros(size1, size2);
            dfdxMatrix = zeros(sdof*NT, sdof*NT);
            for i = 1 : NT
                [fnlt, KT] = AssembleNonlinear_Backbone(U(:,i), eqStructNumeric, i, NT); 
                fnl(:,i) = fnlt;
                index_row_start = sdof*(i-1)+1;
                index_row_end = sdof*i;
                dfdxMatrix(index_row_start:index_row_end, index_row_start:index_row_end) = KT; 
            end
            fnl = reshape(fnl,[],1);
            fnl = FourierInverseTerms.Gamma_p * fnl;
            dfdx = FourierInverseTerms.Gamma_p * dfdxMatrix * FourierInverseTerms.Gamma;
            % ================== AFT end ======================
            R = Om*kron(delta1, B)*Q - fnl - 2*del*kron(I, IalphaFull)*Q;
            dR = [Om*kron(delta1, B)+kron(I, A)-2*del*kron(I, IalphaFull)-dfdx, ...
                kron(delta1, B)*Q, ...
                -2*kron(I, IalphaFull)*Q,...
                Om*kron(delta1, B)*Q + kron(I, A)*Q - 2*del*kron(I, IalphaFull)*Q - dfdx*Q];

    end

    % Scale dynamic force equilibrium (useful for numerical reasons)
    R = 1/fscl*(R);
    dR = 1/fscl*(dR);

    if strcmpi(analysis_type,'nma') || ...
        strcmpi(analysis_type,'nonlinear modal analysis')
        % Scale dynamic force equilibrium by modal amplitude
        % NOTE: We first evaluate the derivative, as we then overwrite R!
        dR(1:end-2,:) = dR(1:end-2,:)/a-R(1:end-2)/a^2*da;
        R(1:end-2) = R(1:end-2)/a;
        
        % Amplitude normalization: The mass of the nonlinear mode shape (all
        % harmonics) is enforced to be one.
        R(end-1) = (Psi'*kron(I,B)*Psi-1);
        dR(end-1,:) = (2*(Psi'*kron(I,B))*dPsi);
        
        % Phase normalization: Velocity of coordinate 'inorm' is enforced to be 
        % zero at t=0.
        R(end) = (1:H)*imag(Psi(inorm+sdof/2+(0:H-1)*sdof));
        dR(end,:) = (1:H)*imag(dPsi(inorm+sdof/2+(0:H-1)*sdof,:));
        % R(end) = (1:H)*imag(Psi(inorm+(sdof:sdof:H*sdof)));
        % dR(end,:) = (1:H)*imag(dPsi(inorm+(sdof:sdof:H*sdof),:));
    end
end