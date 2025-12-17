function B = StrainMatrixNonlinear(G, u)
    du = G*u; 
    A = zeros(6, 9);
    A(1,1:3)=du(1:3,1);
    A(2,4:6)=du(4:6,1);
    A(3,7:9)=du(7:9,1);
    
    A(4,1:3)=du(4:6,1);
    A(4,4:6)=du(1:3,1);
    
    A(5,4:6)=du(7:9,1);
    A(5,7:9)=du(4:6,1);
    
    A(6,1:3)=du(7:9,1);
    A(6,7:9)=du(1:3,1);
    
    B = A * G;

end