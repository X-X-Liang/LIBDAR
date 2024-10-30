function pv = pv_T(T)
% Saturated vapor pressure as a funciton of temperature
% Wagner and Pruss proposed a complex equation for the saturation vapor pressure of water (Wagner and Pruss 1993, 2002). 
% Later this equation formed the basis of the IAPWS formulation.
% Pv is the saturation water vapor pressure (Pa), 
% Tc is the critical point temperature (647.096 K), 
% T is the temperature (K),Array [1 x R] .
% Source: https://journals.ametsoc.org/doi/full/10.1175/JAMC-D-17-0334.1
Tc = 647.096; % K
a = -7.85951783;
b = 1.84408259;
c = -11.7866497;
d = 22.6807411;
e = -15.9618719;
f = 1.80122502;
g = 22064000;
pv = g*exp(Tc./T.*(a*(1-T/Tc)+b*(1-T/Tc).^1.5+c*(1-T/Tc).^3+d*(1-T/Tc).^3.5+e*(1-T/Tc).^4+f*(1-T/Tc).^7.5));
end

