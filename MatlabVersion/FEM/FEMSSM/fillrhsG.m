function Cp = fillrhsG(Mesh,MatParams,Cp,p)
    for p1 = 1 : p-1
        for p2 = 1: p-1
            if (p1+p2) == p
               for i = 1: Cp{p1}.nc
                   for j = 1 : Cp{p2}.nc
                       Avector = Cp{p1}.Avector(i,:) + Cp{p2}.Avector(j,:);
                       pos = find(all(Cp{p}.Avector == repmat(Avector, size(Cp{p}.Avector, 1), 1), 2));
                       if Cp{p}.Corresp(pos)>0
                          vector1 = Cp{p1}.W(:,i);
                          vector2 = Cp{p2}.W(:,j);
                          % Cp{p} = assemble_G(Mesh,MatParams,Cp{p},pos,vector1,vector2);
                          Cp{p} = assemble_G_mex(Mesh,MatParams,Cp{p},pos,vector1,vector2);
                       end
                   end
               end
            end
        end
    end
end