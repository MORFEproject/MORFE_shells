%% ========= Read data from FEM document ===========
clear
clc
restoredefaultpath;
InitModule('RammeShell','SSM');
NodeFile = ".\NumericalExample\ShallowShell_M\NLIST.lis";
ElementFile = ".\NumericalExample\ShallowShell_M\ELIST.lis";
addThicknessColumn(NodeFile, 'TestThickness.txt', 0.01, 0.01);
Mesh = readmesh(NodeFile,ElementFile,'s-s','Average'); 
PlotMesh(Mesh)
NodeFile = 'TestThickness.txt';
Mesh = readmesh(NodeFile,ElementFile,'s-s','Average');
PlotMeshThickness(Mesh,10)
%% ========= Maerial param ==========
E = 70e9;
nu = 0.33;
rho = 2700;
t = Mesh.Node(:,end);
alpha = 0.001*148.9828;
% alpha = 0;
beta = 0;
MatParams = Material(E,nu,t,rho,alpha,beta);
%% ========= Total lagrange pre-calculation ========
PreCalculateParam = PreCalculate(Mesh, MatParams);
Mesh.PreCalculateParam = PreCalculateParam;
%% ========= DPIM param setting ======
MasterMode = [1, 2];
max_order = 3;
max_orderNA = 3;
ComputeMode = 10;
Fmodes = [1];
Fmult = 0.5*[0.1];
Ffreq = 1;
nForce = 1;
omega_mul = 1; 
style = 'c'; % DPIM style (complex norm form:c; real norm form:r; graph style:g)
SSMParams = SSMParam(MasterMode,max_order,max_orderNA,...
    ComputeMode,Fmodes,Fmult,Ffreq,omega_mul,nForce,...
    style,Mesh.FreedofNum);
[K,M] = AssembleKM(Mesh, MatParams);
[Phi, Lambda] = eigs(K(Mesh.Freedof,Mesh.Freedof),M(Mesh.Freedof,Mesh.Freedof),SSMParams.ComputeMode,'sm');
PhiFull = zeros(Mesh.Sdof,ComputeMode);
PhiFull(Mesh.Freedof,:) = Phi;
for i = 1 : SSMParams.ComputeMode
    PhiSelectFull = PhiFull(:, i);
    PhiSelectFull = reshape(PhiSelectFull, 6, []);
    if sum(PhiSelectFull(1:3,:),"all")<0 || abs(min(PhiSelectFull(1:3,:),[],'all'))>abs(max(PhiSelectFull(1:3,:),[],'all'))
       Phi(:,i) = -Phi(:,i);
    end
end
PhiSelectFull = PhiFull(:, 2);
PhiSelectFull = reshape(PhiSelectFull, 6, []);
% PlotDeformation(Mesh, PhiSelectFull', 0.5);
PlotModalShape(Mesh, Phi, ComputeMode, 1, 9); % plot modal shape
Lambda = diag(Lambda);
freq = sqrt(real(Lambda)); 
disp([num2str(freq'./2./pi)])
ForceVector = zeros(Mesh.Sdof,1);
ForceVector(Mesh.ForceDof) = Fmult;
Mesh.ForceVector = ForceVector(Mesh.Freedof); 
%% ========= DPIM routine =====
[Cp, TargetFreq] = ComputesParametrizationHalf(Mesh,SSMParams,MatParams);
Cp = Realification(Cp,SSMParams);
%% ========= Record inverse mapping parameters ====================
howmany = count_terms_dyn(Cp,SSMParams);
[mappings_disp,mappings_vel,mappings_modal_disp,mappings_modal_vel,Avector,fdyn] = ...
    store_dyn_and_map(Mesh,MatParams,SSMParams,Cp,howmany);
%% ========= To calculate the Frequency Response Function using the harmonic balance method ====================
rdyn = write_rdyn(SSMParams,Cp);
% [system, HBParam, Om_HB,Q_HB,Ualpmlitude,~,StableSystem] = ...
%     RunHBMethod(HB_number, Frequency_start, Frequency_end, Arclength_size, Scaling_switch, HB_type);
%     HB_number: Number of HB (Harmonic Balance) superpositions;
%     Frequency_start: FRF calculation starting frequency;
%     Frequency_end: FRF calculation end frequency;
%     Arclength_size: Arc length calculation dimensions;
%     Scaling_switch: Based on the starting frequency for the scaling
%     switch. 0 means off, and 1 means on;
%     HB_type: AdaptiveArclength: Adaptive arc length. 
%              NLvib: See the open-source program by Krack et al. for details.
[system, HBParam, Om_HB,Q_HB,Ualpmlitude,~,StableSystem] = ...
    RunHBMethod(20, TargetFreq(1) * 0.95, TargetFreq(1) * 1.05, 1e-5, 1, 'AdaptiveArclength');
% ref: Krack M, Gross J. Harmonic balance for nonlinear vibration
% problems[M]. Cham: Springer International Publishing, 2019. 
%% ROM's FRF figure
figure
PlotBifurcation(Om_HB/TargetFreq(1), Ualpmlitude,StableSystem, [1], 0.98, 1.02);
%% Inverse mapping figure
U = zeros(size(Ualpmlitude,1)+2, size(Ualpmlitude,2));
U(1:end-2,:) = Ualpmlitude;
U(end-1,:) = 1;
U(end,:) = 1;
[phys_disp1, phys_vel1] = InverseMappingVectorization(U, mappings_disp, mappings_vel, Avector);
[modal_disp1, modal_vel1] = InverseMappingModalVectorization(U, mappings_modal_disp, mappings_modal_vel, Avector);
figure
PlotBifurcation(Om_HB/TargetFreq(1), abs(phys_disp1(74,3,:))/min(t),StableSystem, 1, 0.98, 1.02);
figure
PlotBifurcation(Om_HB/TargetFreq(1), abs(modal_disp1(1,:))/min(t),StableSystem, 1, 0.98, 1.02);
figure
PlotBifurcation(Om_HB/TargetFreq(1), abs(modal_disp1(2,:))/min(t),StableSystem, 1, 0.98, 1.02);
figure
PlotBifurcation(Om_HB/TargetFreq(1), abs(modal_disp1(5,:))/min(t),StableSystem, 1, 0.98, 1.02);
figure
PlotBifurcation(Om_HB/TargetFreq(1), abs(modal_disp1(8,:))/min(t),StableSystem, 1, 0.98, 1.02);



