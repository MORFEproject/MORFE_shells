function HBParam = HBParamSetting
    HBParam.Iter_max = 50000;
%     HBParam.Iter_max = 5e3;
    HBParam.k_max = 15;
    HBParam.cvg_resol = 1e-4;
    HBParam.cvg_begi = 1e0;

    HBParam.loadincr = 3e-5;
    
    HBParam.h_initial = 1e-5;
    HBParam.s_initial = 1;
    HBParam.r_limt = 1.5;
    HBParam.h_limt = 1e-2;    
    HBParam.N_star = 5;    
%     HBParam.freq_start = 0.1;
    HBParam.freq_start = 129.8597 * 1.05;  % 0.28/2/pi  0.8/2/pi
    HBParam.freq_end   = 129.8597 * 1.02;  % 0.70/2/pi  1.2/2/pi
    HBParam.N_H        = 10;       % Number of harmonics 12
    % HBParam.N_T        = 43;
    HBParam.N_T        = 4*HBParam.N_H + 1;      % Number of time samples for IFFT 2^8

end