function Cp = assemble_H(Mesh, MatParams, Cp, pos, vector1, vector2, vector3)
    mult = 1;
    neq=Mesh.FreedofNum;
    F = zeros(Mesh.Sdof,1);
    U1 = zeros(Mesh.Sdof,1);
    U2 = zeros(Mesh.Sdof,1);
    U3 = zeros(Mesh.Sdof,1);
    nhalf = length(vector1)/2;
    U1(Mesh.Freedof) = vector1(1:nhalf);
    U2(Mesh.Freedof) = vector2(1:nhalf);
    U3(Mesh.Freedof) = vector3(1:nhalf);

    ElementNum = Mesh.ElementNum;
    ElementDof = Mesh.ElementDof;
    PreCalculateParam = Mesh.PreCalculateParam;
    for i = 1 : ElementNum
        edof = ElementDof(i,:);
        u1 = U1(edof);
        u2 = U2(edof);
        u3 = U3(edof);

        Fe = Integrate_H(PreCalculateParam, u1, u2, u3, MatParams, i);
        F(edof) = F(edof) + Fe*mult;     
    end
    Cp.rhs(neq+[1:Mesh.FreedofNum],pos) = Cp.rhs(neq+[1:Mesh.FreedofNum],pos) - F(Mesh.Freedof);
end