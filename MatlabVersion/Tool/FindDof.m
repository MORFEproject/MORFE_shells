function dof = FindDof(sdof,freedof)
    freedof_list = [1 : length(freedof)];
    dof = -1*ones(sdof,1);
    dof(freedof) = freedof_list;
end