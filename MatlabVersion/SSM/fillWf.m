function Cp = fillWf(Cp,p,SSMParams)
    for p1 = 2 : p-1
        for p2 = 2 : p-1
            if (p1+p2)==p+1
                for i = 1:Cp{p1}.nc
                    A1=Cp{p1}.Avector(i,:);
                    for j = 1:Cp{p2}.nc
                        A2=Cp{p2}.Avector(j,:);
                        for s = 1:SSMParams.nrom
                           if A1(s)>0
                               Avector=A1+A2;
                               Avector(s)=Avector(s) - 1;
                               pos = find(all(Cp{p}.Avector == repmat(Avector, size(Cp{p}.Avector, 1), 1), 2));
                               Cp{p}.Wf(1:SSMParams.nK,pos)= Cp{p}.Wf(1:SSMParams.nK,pos) + ...
                                   A1(s)*Cp{p1}.W(1:SSMParams.nK,i)*Cp{p2}.f(s,j);
                               Cp{p}.Wf(SSMParams.nK+1:SSMParams.nA,pos) = Cp{p}.Wf(SSMParams.nK+1:SSMParams.nA,pos) + ...
                                   A1(s)*Cp{p1}.W(SSMParams.nK+1:SSMParams.nA,i)*Cp{p2}.f(s,j);
                           end
                        end
                    end
                end
            end
        end
    end
end