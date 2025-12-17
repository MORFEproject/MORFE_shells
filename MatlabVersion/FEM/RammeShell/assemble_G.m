function Cp = assemble_G(Mesh,MatParams,Cp,pos,vector1,vector2)
    mult = 1;
    neq=Mesh.FreedofNum;
    F = zeros(Mesh.Sdof,1);
    U1 = zeros(Mesh.Sdof,1);
    U2 = zeros(Mesh.Sdof,1);
    nhalf = length(vector1)/2;
    U1(Mesh.Freedof) = vector1(1:nhalf);
    U2(Mesh.Freedof) = vector2(1:nhalf);

    ElementNum = Mesh.ElementNum;
    ElementDof = Mesh.ElementDof;
    PreCalculateParam = Mesh.PreCalculateParam;
    for i = 1 : ElementNum
        edof = ElementDof(i,:);
        u1 = U1(edof);
        u2 = U2(edof);

        Fe = Integrate_G(PreCalculateParam, u1, u2, MatParams, i);
        F(edof) = F(edof) + Fe*mult;        
    end
    Cp.rhs(neq+[1:Mesh.FreedofNum],pos) = Cp.rhs(neq+[1:Mesh.FreedofNum],pos) - F(Mesh.Freedof);
end