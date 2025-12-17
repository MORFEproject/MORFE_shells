function Cp = Realification(Cp,SSMParams)
    tic
    disp(['Begin realification process'])
    for p = 1 : SSMParams.max_order
        for i = 1 : Cp{p}.nc
            Avec=Cp{p}.Avector(i,:);
            Ivec=zeros(1,p);
            counter=0;
            for j = 1 : SSMParams.nrom 
                Ivec(counter+1:counter+Avec(j)) = j;
                counter = counter + Avec(j);
            end
            pos = 1;
            coeff = 1.0 + (0.0) * 1i;
            Avec(:) = 0;
            Cp{p} = recursive_C2R(Ivec,p,pos,i,Avec,coeff,Cp{p},SSMParams);
        end
    end
    varInfo = whos('Cp');
    varMemoryMB = varInfo.bytes / 2^20;
    fprintf('Variable "Cp" occupies %.2f MB of memory.\n', varMemoryMB);  
    disp(['End realification process!'])
    toc
    disp(["============================================="])
end