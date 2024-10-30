function sigma  = sigma_T( T )
% Surface tension of water as a function of temperature in Kelvin
% https://en.wikipedia.org/wiki/Surface_tension#cite_ref-IAPWS_31-0
% IAPWS International Association for the Properties of Water and Steam (June 2014). "Revised Release on Surface Tension of Ordinary Water Substance".
% Valid from 273 Kelven to critical point

TC=647.096; % in K
sigma = 235.8*(1-T/TC)^1.256*(1-0.625*(1-T/TC)); % mN/m
sigma = sigma/1000; % N/m

end

