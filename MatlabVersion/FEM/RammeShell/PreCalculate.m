function PreCalculateParam = PreCalculate(Mesh, MatParams)
    [Gauss, Poids, nPtGaus] = GassPoint('2*2*2');
    det_a = zeros(Mesh.ElementNum, nPtGaus); % jacobian matrix
    D = cell(Mesh.ElementNum, nPtGaus); % elastic matrix cell
    R = cell(Mesh.ElementNum, nPtGaus); % covariant matrix cell
    B0 = cell(Mesh.ElementNum, nPtGaus); % linear strain matrix cell
    G = cell(Mesh.ElementNum, nPtGaus); % G matrix cell
    Ba = cell(Mesh.ElementNum, nPtGaus); % EAS strain matrix cell

    invKaa = cell(Mesh.ElementNum, 1); % inverse EAS stiffness matrix cell
    K0 = cell(Mesh.ElementNum, 1); % linear stiffness
    K0a = cell(Mesh.ElementNum, 1); % 
    Ka0 = cell(Mesh.ElementNum, 1); % 
    % linear stiffness (gass points 2*2*2)
    for i = 1 : Mesh.ElementNum
        ElementConnection = Mesh.Element(i,:);
        X = Mesh.Node(ElementConnection,1);
        Y = Mesh.Node(ElementConnection,2);
        Z = Mesh.Node(ElementConnection,3);
        e1 = Mesh.ElementDirection{i,1}; e2 = Mesh.ElementDirection{i,2}; e3 = Mesh.ElementDirection{i,3}; e4 = Mesh.ElementDirection{i,4};
        e5 = Mesh.ElementDirection{i,5}; e6 = Mesh.ElementDirection{i,6}; e7 = Mesh.ElementDirection{i,7}; e8 = Mesh.ElementDirection{i,8};
        VectorDirection = [e1; e2; e3; e4; e5; e6; e7; e8];
        t = MatParams.t(ElementConnection);
        kaa = zeros(4, 4);
        k0 = zeros(48, 48);
        k0a = zeros(48, 4);
        ka0 = zeros(4, 48);
        for j = 1 :  nPtGaus 
            theta1 = Gauss(1, j);
            theta2 = Gauss(2, j); 
            theta3 = Gauss(3, j); 
            N = ShapeFunction(theta1, theta2);
            dN = FirstLocal(theta1,theta2);
            a = LocalBase(dN, N, theta3, X, Y, Z, VectorDirection, t);
            det_a(i, j) = det(a) * Poids(j);
            D{i, j} = ElasticMatrix(MatParams.E, MatParams.nu, a);
            R{i, j} = CovarianteMatrix(a);
            G{i, j} = GradientDeformation(theta3, N, dN, t);
            Ba{i, j} = StrainMatrixAlpha(theta1, theta2, theta3, t, N); 
            B0{i, j} = R{i, j} * G{i ,j};
            kaa = kaa + (Ba{i, j}' * D{i, j} * Ba{i, j} * det_a(i, j));
            k0 = k0 + B0{i, j}' * D{i, j} * B0{i, j} * det_a(i, j);
            k0a = k0a + (B0{i, j}' * D{i, j} * Ba{i, j} * det_a(i, j));
            ka0 = ka0 + (Ba{i, j}' * D{i, j} * B0{i, j} * det_a(i, j));
        end
        invKaa{i, 1} = inv(kaa);
        K0{i, 1} = k0;
        K0a{i, 1} = k0a;
        Ka0{i, 1} = ka0;
    end 
    % mass matrix (gass points 3*3*2)
    [Gauss, Poids, nPtGaus] = GassPoint('3*3*2');
    M = cell(Mesh.ElementNum, 1); 
    for i = 1 : Mesh.ElementNum
        ElementConnection = Mesh.Element(i,:);
        X = Mesh.Node(ElementConnection,1);
        Y = Mesh.Node(ElementConnection,2);
        Z = Mesh.Node(ElementConnection,3);
        e1 = Mesh.ElementDirection{i,1}; e2 = Mesh.ElementDirection{i,2}; e3 = Mesh.ElementDirection{i,3}; e4 = Mesh.ElementDirection{i,4};
        e5 = Mesh.ElementDirection{i,5}; e6 = Mesh.ElementDirection{i,6}; e7 = Mesh.ElementDirection{i,7}; e8 = Mesh.ElementDirection{i,8};
        VectorDirection = [e1; e2; e3; e4; e5; e6; e7; e8];
        t = MatParams.t(ElementConnection);
        m = zeros(48, 48);
        for j = 1 :  nPtGaus 
            theta1 = Gauss(1, j); 
            theta2 = Gauss(2, j); 
            theta3 = Gauss(3, j); 
            N = ShapeFunction(theta1, theta2);
            dN = FirstLocal(theta1,theta2);
            Nmat = ShapeFunctionMatrix(N, theta3, t);
            a = LocalBase(dN, N, theta3, X, Y, Z, VectorDirection, t);
            m = m + MatParams.rho * (Nmat' * Nmat) * det(a) * Poids(j); 
        end
        M{i, 1} = m;
    end 
    PreCalculateParam.det_a = det_a;
    PreCalculateParam.D = D;
    PreCalculateParam.R = R;
    PreCalculateParam.G = G;
    PreCalculateParam.B0 = B0;
    PreCalculateParam.Ba = Ba;
    PreCalculateParam.K0 = K0;
    PreCalculateParam.K0a = K0a;
    PreCalculateParam.Ka0 = Ka0;
    PreCalculateParam.M = M;
    PreCalculateParam.invKaa = invKaa;
end