function Cp = fillrhsH(Mesh,MatParams,Cp,p)
    for p1 = 1 : p-2
        for p2 = 1 : p-2
            for p3 = 1 : p-2
                if (p1+p2+p3)==p
                    for i = 1 : Cp{p1}.nc
                        for j = 1 : Cp{p2}.nc
                            for k = 1 : Cp{p3}.nc
                                Avector=Cp{p1}.Avector(i,:)+Cp{p2}.Avector(j,:)+Cp{p3}.Avector(k,:); 
                                pos = find(all(Cp{p}.Avector == repmat(Avector, size(Cp{p}.Avector, 1), 1), 2));
                                if Cp{p}.Corresp(pos)>0
                                    vector1 = Cp{p1}.W(:,i);
                                    vector2 = Cp{p2}.W(:,j);
                                    vector3 = Cp{p3}.W(:,k);
                                    % Cp{p} = assemble_H(Mesh,MatParams,Cp{p},pos,vector1,vector2,vector3);
                                    Cp{p} = assemble_H_mex(Mesh,MatParams,Cp{p},pos,vector1,vector2,vector3);
                                end
                            end
                        end
                    end
                end
            end
        end
    end
end