function [rho0,c0,varargout] = TaitEOS(p0,varargin)
%% Parse optional input arguments
p = inputParser;
defaultP     = 1e9;    % default value for pressure at bubble wall 1 GPa
addParameter(p,'P',defaultP);
parse(p,varargin{:});   % Parse param.-value pairs
param = p.Results;      % transfer res. to structure
clear p                 % delete the object p that is not needed
P     = param.P;        % Pressure at bubble wall
%% constants
p_STP   = 1e5;       % pressure at standard temperature and pressure (STP)
rho_STP = 998;       % density of water 998kg/m3 at STP
B       = 3.14e8;    % parameter in Tait EOS
n       = 7;         % default 7
c_STP   = 1483;      % speed of sound at STP

% mass density of water at ambient pressure p0 Eq. (3.4) in Ref. LiangJFM2022
rho0    = rho_STP*((p0+B)/(p_STP+B))^(1/n);

% speed of sound in water c0 at ambient pressure p0 Eq. (3.32) in Ref. LiangJFM2022
c0      = c_STP*((p0+B)/(p_STP+B))^((n-1)/(2*n));

% Enthalpy difference between the bubble wall and infinite
H = n/(n-1)*(B+p0)^(1/n)/rho0*((P+B)^((n-1)/n)-(p0+B)^((n-1)/n));  % liquid enthalpy Eq. (3.6) in Liang2022JFM

% sound velocity at the bubble wall
C = (c0^2+(n-1)*H)^0.5;

varargout{1} = H;
varargout{2} = C;

end