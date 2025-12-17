function [M, C, K,delta2,delta1,I] = HBExtensiveMatrix(M, C, K, NH)
    I = eye(2*NH+1,2*NH+1);
    delta1 = zeros(2*NH+1,2*NH+1);
    delta2 = zeros(2*NH+1,2*NH+1);
    for i = 1 : NH
        delta1(2*i, 2*i+1) = i;
        delta1(2*i+1, 2*i) = -i;
        delta2(2*i, 2*i) = -i^2;
        delta2(2*i+1, 2*i+1) = -i^2;
    end
    M = kron(delta2, M);
    C = kron(delta1, C);
    K = kron(I, K);
    M = sparse(M);
    C = sparse(C);
    K = sparse(K);
end