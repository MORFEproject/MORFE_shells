function B_alpha = StrainMatrixAlpha(teta1, teta2, teta3, h, N)
    % thickness interpolation    
    epai=0;
    for k=1:8
        epai=epai+h(k)*N(k);
    end
    
    B_alpha=zeros(6,4);
    
    B_alpha(3,1)=epai/2*teta3;
    B_alpha(3,2)=epai/2*teta3*teta1;
    B_alpha(3,3)=epai/2*teta3*teta2;
    B_alpha(3,4)=epai/2*teta3*teta1*teta2;
end