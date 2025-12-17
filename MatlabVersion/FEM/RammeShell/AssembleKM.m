function [K,M] = AssembleKM(Mesh, MatParams)
    K = sparse(Mesh.Sdof, Mesh.Sdof);
    M = sparse(Mesh.Sdof, Mesh.Sdof);
    ne = Mesh.ElementNum;
    for i = 1 : ne
        edof = Mesh.ElementDof(i, :);
        k = Mesh.PreCalculateParam.K0{i, 1} - ...
            Mesh.PreCalculateParam.K0a{i, 1} * ...
            Mesh.PreCalculateParam.invKaa{i, 1} * ...
            Mesh.PreCalculateParam.Ka0{i, 1};
        m = Mesh.PreCalculateParam.M{i, 1};
        K(edof, edof) = K(edof, edof) + k;
        M(edof, edof) = M(edof, edof) + m;
    end
end