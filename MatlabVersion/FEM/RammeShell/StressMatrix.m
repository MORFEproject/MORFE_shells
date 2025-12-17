function M = StressMatrix(SIG_M)
    M=zeros(9,9);
    for i=1:3
         M(i,i)=SIG_M(1);
         M(i,i+3)=SIG_M(4);
         M(i,i+6)=SIG_M(6);
         M(i+3,i)=SIG_M(4);
         M(i+3,i+3)=SIG_M(2);
         M(i+3,i+6)=SIG_M(5);
         M(i+6,i)=SIG_M(6);
         M(i+6,i+3)=SIG_M(5);
         M(i+6,i+6)=SIG_M(3);
    end
end