function G = GradientDeformation(teta3, N, dN, h)
    G=zeros(9,6*8);   
    for  K=1:8
        indice=6*(K-1);
        TET3=teta3*0.5*h(K);
        G(1,1+indice)=dN(K,1);
        G(1,4+indice)=TET3*dN(K,1);
        G(2,2+indice)=G(1,1+indice);
        G(2,5+indice)=G(1,4+indice);
        G(3,3+indice)=G(1,1+indice);
        G(3,6+indice)=G(1,4+indice);
        G(4,1+indice)=dN(K,2);
        G(4,4+indice)=TET3*dN(K,2);
        G(5,2+indice)=G(4,1+indice);
        G(5,5+indice)=G(4,4+indice);
        G(6,3+indice)=G(4,1+indice);
        G(6,6+indice)=G(4,4+indice);
        G(7,4+indice)=0.5*h(K)*N(K);
        G(8,5+indice)=G(7,4+indice);
        G(9,6+indice)=G(7,4+indice);    
    end    
end