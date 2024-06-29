clear all, clc
beta = 7;
%%  Data
V0    = (50-beta)/5;      % [L], volume
M0    = 0.005;     % [mol/L], initial concentration of monomer
S0    = 0.1;       % [mol/L], initial concentration of solvent
I0    = 0.02;    % [mol/L]
X0    = 0;
MW_M  = 120;     % [g/mol], molecular weight monomer
MW_S  = 85;      % [g/mol], molecular weight solvent
MW_I  = 270;     % [g/mol], molecular weight initiator
f     = 0.5;     % [-]
kd    = 4e-6; % [s-1]
kp0   = 750;     % [L/mol/s]
kt0   = 9e6;   % [L/mol/s]
kpD0  = 3.7e11;    % [L/mol/s]
ktD0  = 2.8e8;     % [L/mol/s]
Cn    = 25;      % [-]
CRD   = 180;     % [-]
X0 = 0;          % [-] Initial Conversion

%% Resolution
mM0 = M0*V0*MW_M;
mS0 = S0*V0*MW_S;
mI0 = I0*MW_I;

Cfm = 0.01;
Cfs = 0.001;
Ct = 1000;

tspan = [0:0.1:2*3600]';
C0 = [I0 M0 S0 X0]';
[t,C] = ode23s(@(t,C)Batch_Diffusion(t,C,f,kd,V0,M0,MW_M,mM0,mS0,MW_S,kp0,kt0,kpD0,ktD0,Cn,Cfm,Cfs,CRD),tspan,C0);
I = C(:,1);
M = C(:,2);
S = C(:,3);
X = C(:,4);

% Mass Calculations
mS = S.*MW_S.*V0;
mM = M0.*(1-X).*MW_M.*V0;
mP = mM0-mM +mS0 - mS;
wP = mP./(mP+mS+mM);

% Kinetic Constants
kp = (1/kp0+exp(Cn.*wP)./kpD0).^(-1);
kt = (1/kt0+exp(Cn.*wP)./ktD0).^(-1)+CRD*kp.*(1-wP);
kfm = kp*Cfm;
kfs = kp*Cfs;
ktc = kt./(1+Ct);
ktd = Ct.*kt./(1+Ct);

R = sqrt((2.*f.*kd.*I)./kt);

% Characteristic Times
taup = 1./(kp.*M);
tautc = 1./(ktc.*R);
tautd = 1./(ktd.*R);
taufm = 1./(kfm.*M);
taufs = 1./(kfs.*S);

%Evaluate alpha
beta = taup./tautc;
gamma = taup./taufm+taup./tautd+taup./taufs;
alpha = beta + gamma;

% Instantaneous Properties
mu_1 = 1./(alpha+beta./2);
mu_2 = (2.*gamma+3.*beta)./(alpha.^2.*(gamma+0.5*beta));

% Polymer Instantaneous Properties
Mn_inst_av = MW_M.*mu_1;
Mw_inst_av = MW_M.*mu_2./mu_1;

disp(Mw_inst_av(end))

%% Function 1 - With Diffusion
function dF = Batch_Diffusion(t,C,f,kd,V0,M0,MW_M,mM0,mS0,MW_S,kp0,kt0,kpD0,ktD0,Cn,Cfm,Cfs,CRD)
% Unknowns

I = C(1);
M = C(2);
S = C(3);
X = C(4);

% Preliminar Calculations
% Calculations
mM   = M0.*(1-X).*MW_M.*V0;     %[g], concentrazione*volume*peso molecolare=massa, concentrazione moltiplicata per (1-X) perchè non valutiamo le condizioni iniziali                
mS   = S.*MW_S.*V0;     %[g]
mP   = mM0-mM +mS0 - mS;                  %[g], il polimero viene prodotto e all'inizio è 0 quindi la sua massa sarà la differenza tra la quantità di monomero iniziale e attuale
wP   = mP./(mP+mS+mM);          %[-], frazione massiva NO mI

% Kinetic Constants
kp    = (1/kp0+exp(Cn.*wP)./kpD0).^(-1);                  %[L/mol/s] formula nel testo
kt    = ((1/kt0+exp(Cn.*wP)./ktD0).^(-1))+CRD*kp.*(1-wP); %[L/mol/s] formula nel testo
kfm   = kp.*Cfm;                                          %[L/mol/s] formula nel testo
kfs   = kp.*Cfs;                                          %[L/mol/s]

R = sqrt((2.*f.*kd.*I)./kt);

dI = -kd.*I;
dM = -(kp+kfm).*M.*R;
dS = -kfs.*S.*R;
dX = (kp+kfm).*(1-X).*R;

dF = [dI dM dS dX]';
end
