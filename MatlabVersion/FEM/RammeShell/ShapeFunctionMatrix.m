function N = ShapeFunctionMatrix(Ncomponents, theta3, h)
    Nmat1 = SubMatrix(Ncomponents(1),theta3,h(1));
    Nmat2 = SubMatrix(Ncomponents(2),theta3,h(2));
    Nmat3 = SubMatrix(Ncomponents(3),theta3,h(3));
    Nmat4 = SubMatrix(Ncomponents(4),theta3,h(4));

    Nmat5 = SubMatrix(Ncomponents(5),theta3,h(5));
    Nmat6 = SubMatrix(Ncomponents(6),theta3,h(6));
    Nmat7 = SubMatrix(Ncomponents(7),theta3,h(7));
    Nmat8 = SubMatrix(Ncomponents(8),theta3,h(8));

    N = [Nmat1, Nmat2, Nmat3, Nmat4, Nmat5, Nmat6, Nmat7, Nmat8];

    function Nmat = SubMatrix(N, theta3, h)
        Nmat = zeros(3, 6);
        Nmat(1, 1) = N; Nmat(2, 2) = N; Nmat(3, 3) = N;
        Nmat(1, 4) = 0.5*h*theta3*N; Nmat(2, 5) = 0.5*h*theta3*N; Nmat(3, 6) = 0.5*h*theta3*N;
    end
end