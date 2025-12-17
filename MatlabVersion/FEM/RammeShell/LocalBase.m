function a = LocalBase(dN, N, theta3, X, Y, Z, VectorDirection, h)
    a=zeros(3,3);
    for k = 1 : 8 
        a(1,1)=a(1,1)+dN(k,1)*X(k)+0.5*theta3*h(k)*dN(k,1)*VectorDirection(k,1);
        a(1,2)=a(1,2)+dN(k,1)*Y(k)+0.5*theta3*h(k)*dN(k,1)*VectorDirection(k,2);
        a(1,3)=a(1,3)+dN(k,1)*Z(k)+0.5*theta3*h(k)*dN(k,1)*VectorDirection(k,3);
        
        a(2,1)=a(2,1)+dN(k,2)*X(k)+0.5*theta3*h(k)*dN(k,2)*VectorDirection(k,1);
        a(2,2)=a(2,2)+dN(k,2)*Y(k)+0.5*theta3*h(k)*dN(k,2)*VectorDirection(k,2);
        a(2,3)=a(2,3)+dN(k,2)*Z(k)+0.5*theta3*h(k)*dN(k,2)*VectorDirection(k,3);
        
        a(3,1)=a(3,1)+0.5*h(k)*N(k)*VectorDirection(k,1);
        a(3,2)=a(3,2)+0.5*h(k)*N(k)*VectorDirection(k,2);
        a(3,3)=a(3,3)+0.5*h(k)*N(k)*VectorDirection(k,3);        
    end
end