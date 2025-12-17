function FourierInverse = DPIMFourierInverse(sdof, HBParam)
    NT = HBParam.N_T; 
    NH = HBParam.N_H; 
    I = eye(sdof,sdof);

    Ece = sparse(NT, 2*NH+1);
    Ece(:, 1) = 1;
    for i = 1 : NT
        for j = 1 : NH
            Ece(i, 2*j) = cos(2*pi*(i-1)*j/NT);
            Ece(i, 2*j+1) = sin(2*pi*(i-1)*j/NT);
        end
    end
    Gamma = kron(Ece, I);
    Gamma = sparse(Gamma);

    Ece_inv = sparse(2*NH+1, NT);
    Ece_inv(1, :) = 0.5;
    for i = 1 : NH
        for j = 1 : NT
            Ece_inv(2*i, j) = cos(2*pi*i*(j-1)/NT);
            Ece_inv(2*i+1, j) = sin(2*pi*i*(j-1)/NT);
        end
    end
    Ece_inv = Ece_inv * (2/NT);
    Gamma_p = kron(Ece_inv, I);
    Gamma_p = sparse(Gamma_p);

    FourierInverse.Gamma = Gamma;
    FourierInverse.Gamma_p = Gamma_p;
end