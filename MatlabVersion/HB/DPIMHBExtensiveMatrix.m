function [delta1, I] = DPIMHBExtensiveMatrix(NH)
    I = eye(2*NH+1,2*NH+1);
    delta1 = zeros(2*NH+1,2*NH+1);
    for i = 1 : NH
        delta1(2*i, 2*i+1) = i;
        delta1(2*i+1, 2*i) = -i;
    end
end