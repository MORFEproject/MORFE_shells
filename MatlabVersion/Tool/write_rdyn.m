function rdyn1 = write_rdyn(SSMParams,Cp)
  tic
  disp(['Begin write Reduced dynamic equations!'])
  disp(['norm form = ', num2str(SSMParams.nrom)])
  rdyn = cell(SSMParams.nz,1);
  rdyn1 = cell(SSMParams.nz,1);
  for i = 1 : SSMParams.nz
      str1 = "z";
      str2 = num2str(i);
      str3 = "' =";
      str = str1+str2+str3;
      rdyn{i} = str;
      str4 = " ";
      rdyn1{i} = str4;
  end
  for p = 1:SSMParams.max_order
      for c = 1:Cp{p}.nc
          Avector=Cp{p}.Avector(c,:);
          monomial = "";
          for d = 1:SSMParams.nrom 
              if (Avector(d)~=0)
                  monomial = monomial+"*z"+num2str(d)+"^"+num2str(Avector(d));
              end
          end
          for j = 1 : SSMParams.nm
              rcoeff=2*real(Cp{p}.fr(j,c));
              icoeff=-2*imag(Cp{p}.fr(j,c));
              if abs(rcoeff)>1e-20
                  rdyn{j} =  rdyn{j} + " + " + num2str(rcoeff) + monomial;
                  rdyn1{j} = rdyn1{j} + " + " + num2str(rcoeff) + monomial;
              end
              if abs(icoeff)>1e-20
                  rdyn{j+SSMParams.nm} = ...
                      rdyn{j+SSMParams.nm} + " + " + num2str(icoeff) + monomial;
                  rdyn1{j+SSMParams.nm} = ...
                      rdyn1{j+SSMParams.nm} + " + " + num2str(icoeff) + monomial;
              end
          end
      end
  end
  disp(['Reduced dynamic equations:'])
  % write in text
  folderName = 'Output';
  currentFolder = pwd;
  targetFolderFullPath = fullfile(currentFolder, folderName);
  pathString = genpath(targetFolderFullPath);
  addpath(pathString);
  filePath = 'Output\Equations.txt';
  fileID = fopen(filePath, 'w');
  % 
  for i = 1:length(rdyn)
      fprintf(fileID, '%s\n', rdyn{i});
  end
  % 
  fclose(fileID);
  % Processing the ROM (Reduced-Order Model) equations for subsequent HB (Harmonic Balance) calculations
  EquationSeparte_numeric(filePath)
  toc
  disp(['End write Reduced dynamic equations!'])
  disp(["============================================="])
end