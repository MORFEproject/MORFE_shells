function [KL, KT] = AssembleNonlinear(Ureduce, Mesh, MatParams)
    U = zeros(Mesh.Sdof,1);
    U(Mesh.Freedof) = Ureduce;
    element = Mesh.Element;
    node = Mesh.Node;
    KT = sparse(Mesh.Sdof, Mesh.Sdof);
    KL = sparse(Mesh.Sdof, Mesh.Sdof);   
    for i = 1 : Mesh.ElementNum
        edof = Mesh.ElementDof(i, :);
        u = U(edof);
        [kT, kNL] = ElementStiffnessNonlinear(i, u, Mesh);
        KL(edof, edof) = KL(edof, edof) + kNL;
        KT(edof, edof) = KT(edof, edof) + kT;
    end
end