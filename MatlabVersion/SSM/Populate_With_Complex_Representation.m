function Cp = Populate_With_Complex_Representation(Cp, SSMParams)
    tic;
    disp('Begin populating coefficients with complex representation (identity mapping)');

    for p = 1 : SSMParams.max_order
        Cp{p}.Wr = Cp{p}.W;
        Cp{p}.fr = Cp{p}.f;
    end
    
    varInfo = whos('Cp');
    varMemoryMB = varInfo.bytes / 2^20;
    fprintf('Variable "Cp" occupies %.2f MB of memory.\n', varMemoryMB);
    disp('End populating coefficients!');
    toc;
    disp("=============================================");
end