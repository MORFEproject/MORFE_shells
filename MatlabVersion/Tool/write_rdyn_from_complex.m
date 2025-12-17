function rdyn_z = write_rdyn_from_complex(SSMParams, Cp)
    tic;
    disp('Begin writing complex reduced dynamic equations!');
    
    nm = SSMParams.nm;
    
    rdyn_z = cell(nm, 1);
    for i = 1 : nm
        rdyn_z{i} = "z" + num2str(i) + "' =";
    end

    % 
    for p = 1 : SSMParams.max_order
        for c = 1 : Cp{p}.nc
            Avector = Cp{p}.Avector(c, :);
            monomial_str = "";
            for d = 1 : SSMParams.nrom 
                if (Avector(d) ~= 0)
                    if d <= nm
                        var_name = "z" + num2str(d);
                    else
                        var_name = "z_bar" + num2str(d - nm);
                    end
                    
                    monomial_str = monomial_str + "*" + var_name;
                    if Avector(d) > 1
                        monomial_str = monomial_str + "^" + num2str(Avector(d));
                    end
                end
            end

            for j = 1 : nm
                coeff = Cp{p}.fr(j, c);
                if abs(coeff) > 1e-20
                    coeff_str = sprintf('(%f + %fi)', real(coeff), imag(coeff)); 
                    rdyn_z{j} = rdyn_z{j} + " + " + coeff_str + monomial_str;
                end
            end
        end
    end

    folderName = 'Output';
    if ~exist(folderName, 'dir')
       mkdir(folderName)
    end
    filePath = fullfile(folderName, 'Equations_complex.txt');
    
    fileID = fopen(filePath, 'w');
    for i = 1:length(rdyn_z)
        fprintf(fileID, '%s\n', rdyn_z{i});
    end
    fclose(fileID);

    toc;
    disp(['End writing complex reduced dynamic equations to ' filePath]);
    disp("=============================================");
end