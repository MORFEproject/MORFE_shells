function Cp = InitParametrization(SSMParams)
    Cp = cell(SSMParams.max_order+1,1);
    for i = 1 : SSMParams.max_order
        [Avector, Corresp, nc, ncindep] = indexset(SSMParams, i);
        Cp{i}.Avector = Avector;
        Cp{i}.Corresp = Corresp;
        Cp{i}.nc = nc;
        Cp{i}.ncindep = ncindep;
        Cp{i}.W = zeros(SSMParams.nA,Cp{i}.nc);
        Cp{i}.Wr = zeros(SSMParams.nA,Cp{i}.nc);
        Cp{i}.f = zeros(SSMParams.nrom,Cp{i}.nc);
        Cp{i}.Wf = zeros(SSMParams.nA,Cp{i}.nc);
        Cp{i}.fr = zeros(SSMParams.nrom,Cp{i}.nc);
        Cp{i}.rhs = zeros(SSMParams.nA,Cp{i}.nc);
    end
    
    function [a, corresp, lena, nc] = indexset(SSMParams, i)
        nz = SSMParams.nz; 
        nzf = 2 * SSMParams.nForce;
        nzhalf = nz/2;
        nzfhalf = nzf/2;
        ndof = SSMParams.nrom;
        a = generateMonomialExponents(ndof, i); 
        lena = length(a);
        corresp=zeros(lena,1);
        nc = 0;
        for i = 1 : lena
            acomponent = a(i,:);
            orderna = sum(acomponent(nz+1:ndof));
            if orderna <= SSMParams.max_orderNA
                b = [acomponent(nzhalf+1:nz), acomponent(1:nzhalf), acomponent(nz+nzfhalf+1:ndof),acomponent(nz+1:nz+nzfhalf)];
                if corresp(i) == 0
                    nc = nc + 1;
                    for j = 1 : lena
                       if isequal(a(j,:),b)
                           corresp(j) = -i;
                           corresp(i) = j;
                           continue
                       end
                    end
                end 
            end
        end
    end
end