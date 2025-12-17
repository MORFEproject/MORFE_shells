function SSMParams = SSMParamInit(Mesh)
    MasterMode = [1];
    max_order = 5;
    max_orderNA = 5; 
    ComputeMode = 10;
    Fmodes = [1];
    Fmult = 0.5*[5];
    Ffreq = 1;
    nForce = 1; 
    omega_mul = 1; 
    style = 'g';
    SSMParams = SSMParam(MasterMode,max_order,max_orderNA,ComputeMode,...
        Fmodes,Fmult,Ffreq,omega_mul,nForce,style,Mesh.FreedofNum);
end