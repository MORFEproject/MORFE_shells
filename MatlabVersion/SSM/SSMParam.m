classdef SSMParam
   properties (SetAccess = public)
       MasterMode; 
       max_order; 
       max_orderNA; 
       ComputeMode; 
       Fmodes; 
       Fmult;
       Ffreq; 
       omega_mul;
       nForce; 
       style;
       nK; 
       nA;
       nrom; 
       nz; 
       nMat; 
       nm; 

   end
   methods 
       function obj = SSMParam(MasterMode,max_order,max_orderNA,ComputeMode,Fmodes,Fmult,Ffreq,omega_mul,nForce,style,Freedof)
          obj.MasterMode = MasterMode;
          obj.max_order = max_order;
          obj.max_orderNA = max_orderNA; 
          obj.ComputeMode = ComputeMode;
          obj.Fmodes = Fmodes;
          obj.Fmult = Fmult;
          obj.nForce = nForce;
          obj.omega_mul = omega_mul;
          obj.Ffreq = Ffreq;
          obj.style = style;
          obj = obj.NumFirstSystemAfterMapping();
          obj = obj.NumFirstSystemBeforeMapping(Freedof);
          obj = obj.NumMasterModeFirstSystemAfterMapping();
          obj = obj.NumMasterModeSecondSystemAfterMapping();
          obj = obj.NumTotal();
       end
       function obj = NumFirstSystemBeforeMapping(obj,FreedofNum)
          obj.nK =  FreedofNum;
          obj.nA = 2 * obj.nK;
       end
       function obj = NumFirstSystemAfterMapping(obj)
          obj.nrom = 2*obj.nForce + 2*length(obj.MasterMode);
       end
       function obj = NumMasterModeFirstSystemAfterMapping(obj)
          obj.nz =  2*length(obj.MasterMode);
       end
       function obj = NumMasterModeSecondSystemAfterMapping(obj)
          obj.nm =  length(obj.MasterMode);
       end       
       function obj = NumTotal(obj)
          obj.nMat =  obj.nA + 2*length(obj.MasterMode);
       end
   end
    
end