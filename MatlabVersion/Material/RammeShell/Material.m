classdef Material
   properties (SetAccess = private)
       E;
       nu;
       rho;
       t;
       alpha;
       beta;
   end
   methods
       function obj = Material(E,nu,t,rho,alpha,beta)
           obj.E = E;
           obj.nu = nu;
           obj.t = t;
           obj.rho = rho;
           obj.alpha = alpha;
           obj.beta = beta;
       end
   end
end