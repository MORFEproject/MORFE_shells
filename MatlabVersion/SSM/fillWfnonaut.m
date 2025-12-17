function Cp = fillWfnonaut(Cp,p,Avector,ind_rhs,SSMParams)
    for r = SSMParams.nz+1:SSMParams.nrom
       if  Avector(r)>0
           for s = 1:SSMParams.nz
               fs_r = Cp{1}.f(s,r);
               if abs(fs_r)>10^(-8)    % only fills if the reduced dyn is nzero
                   Av_W=Avector;
                   Av_W(s) = Av_W(s) + 1;
                   Av_W(r) = Av_W(r) - 1;
                   ind_W = find(all(Cp{p}.Avector == repmat(Av_W, size(Cp{p}.Avector, 1), 1), 2));
                   Cp{p}.Wf(1:SSMParams.nK,ind_rhs) = Cp{p}.Wf(1:SSMParams.nK,ind_rhs) +...
                       Cp{p}.W(1:SSMParams.nK,ind_W)*fs_r*Av_W(s);
                   Cp{p}.Wf(SSMParams.nK+1:SSMParams.nA,ind_rhs)= Cp{p}.Wf(SSMParams.nK+1:SSMParams.nA,ind_rhs)+...
                       Cp{p}.W(SSMParams.nK+1:SSMParams.nA,ind_W)*fs_r*Av_W(s);
               end
           end
       end
    end
end
