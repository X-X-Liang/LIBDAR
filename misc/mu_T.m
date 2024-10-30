function mu = mu_T( T )
%  dynamic viscosity of water as a function of temperature in Kelvin
%  https://en.wikipedia.org/wiki/Viscosity#Water
%  Equations Ref: https://www.engineersedge.com/physics/water__density_viscosity_specific_weight_13146.htm
    A=2.414e-5; % Pa.s
    B=247.8;    % K
    C=140;      % K

 mu = A*10^(B/(T-C));
end

