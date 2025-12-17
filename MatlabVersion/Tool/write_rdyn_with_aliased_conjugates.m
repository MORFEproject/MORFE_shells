function rdyn_z = write_rdyn_with_aliased_conjugates(SSMParams, Cp)
    tic;
    disp('Begin writing aliased-conjugate reduced dynamic equations!');
    
    nm = SSMParams.nm;
    nz = SSMParams.nz; % nz = 2*nm
    
    rdyn_z = cell(nz, 1);
    for i = 1 : nz
        rdyn_z{i} = "z" + num2str(i) + "' =";
    end

    for p = 1 : SSMParams.max_order
        for c = 1 : Cp{p}.nc
            Avector = Cp{p}.Avector(c, :);
            for j = 1 : nm
                coeff = Cp{p}.fr(j, c);

                if abs(coeff) < 1e-20
                    continue;
                end

                monomial_str = "";
                for d = 1 : SSMParams.nrom % nrom = nz
                    if (Avector(d) ~= 0)
                        var_name = "z" + num2str(d);
                        monomial_str = monomial_str + "*" + var_name;
                        if Avector(d) > 1
                            monomial_str = monomial_str + "^" + num2str(Avector(d));
                        end
                    end
                end
                coeff_str = sprintf('(%f + %fi)', real(coeff), imag(coeff));
                rdyn_z{j} = rdyn_z{j} + " + " + coeff_str + monomial_str;

                conj_coeff = conj(coeff);
                conj_coeff_str = sprintf('(%f + %fi)', real(conj_coeff), imag(conj_coeff));
                
                monomial_conj_str = "";
                for d = 1 : SSMParams.nrom
                    if (Avector(d) ~= 0)
                        if d <= nm
                            var_name = "z" + num2str(d + nm);
                        else
                            var_name = "z" + num2str(d - nm);
                        end
                        monomial_conj_str = monomial_conj_str + "*" + var_name;
                        if Avector(d) > 1
                            monomial_conj_str = monomial_conj_str + "^" + num2str(Avector(d));
                        end
                    end
                end
                rdyn_z{j+nm} = rdyn_z{j+nm} + " + " + conj_coeff_str + monomial_conj_str;
            end
        end
    end

    folderName = 'Output';
    if ~exist(folderName, 'dir'), mkdir(folderName); end
    filePath = fullfile(folderName, 'Equations_aliased.txt');
    
    fileID = fopen(filePath, 'w');
    for i = 1:length(rdyn_z)
        fprintf(fileID, '%s\n', rdyn_z{i});
    end
    fclose(fileID);

    toc;
    disp(['End writing aliased-conjugate reduced dynamic equations to ' filePath]);
    disp("=============================================");
end