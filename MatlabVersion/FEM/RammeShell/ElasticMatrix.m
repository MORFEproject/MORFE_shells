function D = ElasticMatrix(E, nu, a)
    % Tenseur metrique [aa] de [a], [a] etant base quelconque
    aa=a*a';
    % [inv_aa]: Base contravariante de la base [a]
    aa_inv=inv(aa);
    %
    YG=E;
    NU=nu;
    %
    LAMDA=YG*NU/((1.D0 + NU)*(1.D0-2.D0*NU));
    MU=YG/(2.D0*(1.D0 + NU));
    LAMDA2MU=LAMDA+2.D0*MU;
    
    D=zeros(6,6);
    
    D(1,1)=LAMDA2MU*aa_inv(1,1)*aa_inv(1,1);
    D(1,2)=LAMDA*aa_inv(1,1)*aa_inv(2,2)+2.D0*MU*aa_inv(1,2)*aa_inv(1,2);
    D(1,3)=LAMDA*aa_inv(1,1)*aa_inv(3,3)+2.D0*MU*aa_inv(1,3)*aa_inv(1,3);
    D(1,4)=LAMDA2MU*aa_inv(1,1)*aa_inv(1,2);
    D(1,5)=LAMDA*aa_inv(1,1)*aa_inv(2,3)+2.D0*MU*aa_inv(1,2)*aa_inv(1,3);
    D(1,6)=LAMDA2MU*aa_inv(1,1)*aa_inv(1,3);
    
    D(2,1)=D(1,2);
    D(2,2)=LAMDA2MU*aa_inv(2,2)*aa_inv(2,2);
    D(2,3)=LAMDA*aa_inv(2,2)*aa_inv(3,3)+2.D0*MU*aa_inv(2,3)*aa_inv(2,3);
    D(2,4)=LAMDA2MU*aa_inv(2,2)*aa_inv(2,1);
    D(2,5)=LAMDA2MU*aa_inv(2,2)*aa_inv(2,3);
    D(2,6)=LAMDA*aa_inv(2,2)*aa_inv(3,1)+2.D0*MU*aa_inv(2,3)*aa_inv(2,1);
    
    D(3,1)=D(1,3);
    D(3,2)=D(2,3);
    D(3,3)=LAMDA2MU*aa_inv(3,3)*aa_inv(3,3);
    D(3,4)=LAMDA*aa_inv(3,3)*aa_inv(1,2)+2.D0*MU*aa_inv(3,1)*aa_inv(3,2);
    D(3,5)=LAMDA2MU*aa_inv(3,3)*aa_inv(3,2);
    D(3,6)=LAMDA2MU*aa_inv(3,3)*aa_inv(3,1);
    
    D(4,1)=D(1,4);
    D(4,2)=D(2,4);
    D(4,3)=D(3,4);
    D(4,4)=(LAMDA+MU)*aa_inv(1,2)*aa_inv(1,2)+MU*aa_inv(1,1)*aa_inv(2,2);
    D(4,5)=(LAMDA+MU)*aa_inv(1,2)*aa_inv(2,3)+MU*aa_inv(1,3)*aa_inv(2,2);
    D(4,6)=(LAMDA+MU)*aa_inv(1,2)*aa_inv(3,1)+MU*aa_inv(1,1)*aa_inv(2,3);
    
    D(5,1)=D(1,5);
    D(5,2)=D(2,5);
    D(5,3)=D(3,5);
    D(5,4)=D(4,5);
    D(5,5)=(LAMDA+MU)*aa_inv(2,3)*aa_inv(2,3)+MU*aa_inv(2,2)*aa_inv(3,3);
    D(5,6)=(LAMDA+MU)*aa_inv(2,3)*aa_inv(3,1)+MU*aa_inv(2,1)*aa_inv(3,3);
    
    D(6,1)=D(1,6);
    D(6,2)=D(2,6);
    D(6,3)=D(3,6);
    D(6,4)=D(4,6);
    D(6,5)=D(5,6);
    D(6,6)=(LAMDA+MU)*aa_inv(3,1)*aa_inv(3,1)+MU*aa_inv(3,3)*aa_inv(1,1);    


end