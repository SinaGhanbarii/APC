%% Excercise - 5
beta = 6;
V0 = (50-beta)/5;        % [L]
M0 = 5e-3;      % [mol/lit]
I0 = 0.1;      % [mol/lit]
MW_M = 120;     % [gr/mol]
MW_S = 85;      % [gr/mol]
MW_I = 270;     % [gr/mol]
f = 0.5;        % [-]
kd = 4e-6;   % [s^-1]
kp0 = 700;      % [lit/mol/s]
kt0 = 9e+3;    % [lit/mol/s]
kpD0 = 3.7e+11;   % [lit/mol/s]
ktD0 = 2.8e5;     % [lit/mol/s]
Cn = 25;        % [-]
CRD = 180;      % [-]

%% Resolution
% Preliminar Calculations
% mS0 = (V0-mM0/rho_M)*rho_S; % [gr]
% mI0 = w_I0*mM0;     % [gr]
% I0 = (mI0/MW_I)/V0; % [mol/lit]
mI0 = I0*V0*MW_I;
mM0 = M0*V0*MW_M;

%% Function 1 --- Diffusion Case
function dF = Batch_Diffusion(X,C,f,kd,V0,M0,MW_M,mM0,mS0,kp0,kt0,kpD0,ktD0,Cn,Cfm,CRD)
% Unknowns 
I = C(1);
M = C(2);
% Calculations
mM = M0.*(1-X).*MW_M.*V0;
mP = mM0 - mM;          % [gr]
wP = mP./(mP+mS0+mM);   % [-]
% Kinetics Constants 
kp = (1/kp0+exp(Cn.*wP)./kpD0).^(-1);   % [lit/mol/s]
kt = ((1/kt0+exp(Cn.*wP)./ktD0).^(-1))+CRD*kp.*(1-wP);
kfm = kp.*Cfm;  % [lit/mol/s]

dI = -(kd.*I)./((kp+kfm+kfs).*(1-X).*sqrt(2.*f.*kd.*I./kt));
dM = -M./(1-X);

dF = [dI dM]';
end