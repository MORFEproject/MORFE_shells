function [Cp] = DPIMParamHigherOrderHalf(Mesh,SSMParams,MatParams,Cp,freq,Phi,...
    K,M,C,Sol,Rhs,Mat,Tri,AY,XTA)
    disp(["Higher orders"])
    for p = 2 : SSMParams.max_order
        disp(["============================================="])
        disp(['order = ',num2str(p)])
        Cp = fillrhsG(Mesh,MatParams,Cp,p);
        Cp = fillrhsH(Mesh,MatParams,Cp,p);
        Cp = fillWf(Cp,p,SSMParams);
        for i = 1:Cp{p}.nc
            corresp = Cp{p}.Corresp(i);
            if corresp>0
                Cp = fillWfnonaut(Cp,p,Cp{p}.Avector(i,:),i,SSMParams);
                [Sol,Rhs,Mat,Cp] = ...
                     homological_HALF(Rhs,Mat,Tri,Cp,p,i,M,K,C,AY,SSMParams,MatParams);
            elseif corresp<0 
                Cp{p}.W(:,i)=conj(Cp{p}.W(:,-corresp));
                Cp{p}.f(1:SSMParams.nm,i)=conj(Cp{p}.f(SSMParams.nm+1:2*SSMParams.nm,-corresp));
                Cp{p}.f(SSMParams.nm+1:2*SSMParams.nm,i)=conj(Cp{p}.f(1:SSMParams.nm,-corresp));
            end
        end
    end
end