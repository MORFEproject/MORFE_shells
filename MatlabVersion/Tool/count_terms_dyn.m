function howmany = count_terms_dyn(Cp,SSMParams)
    howmany = 0;
    for p = 1:SSMParams.max_order
       for  i = 1:Cp{p}.nc 
           howmany = howmany + 1; 
       end
    end
end