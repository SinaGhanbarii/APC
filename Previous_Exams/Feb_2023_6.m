clear all, clc
beta = 6;
%%  Data
V0    = (60-beta)/4;      % [L], volume
M0    = 5*(1+beta);     % [mol/L], initial concentration of monomer
S0    = 10;       % [mol/L], initial concentration of solvent
I0    = 0.01;    % [mol/L]
X0    = 0;
MW_M  = 120;     % [g/mol], molecular weight monomer
MW_S  = 85;      % [g/mol], molecular weight solvent
MW_I  = 270;     % [g/mol], molecular weight initiator
% The dinsity values hasn't give in the problem data!
rho_M = 940;     % [g/L], monomer density
rho_S = 880;     % [g/L], solvent density
f     = 0.5;     % [-]
kd    = 4e-6; % [s-1]
kp0   = 700;     % [L/mol/s]
kt0   = 9e6;   % [L/mol/s]
kpD0  = 3.7e11;    % [L/mol/s]
ktD0  = 2.8e8;     % [L/mol/s]
Cn    = 25;      % [-]
CRD   = 180;     % [-]

%% Resolution
mM0 = M0*V0*MW_M;
mS0 = (V0-mM0/rho_M)*rho_S;
mI0 = I0*MW_I;

Cfm = 0.01;
Cfs = 0.003;
Ct = 1000;
C0 = [I0 M0 S0 X0]';
tau = [0:0.1:4*3600]';
%% RESOLUTION WITH DIFFUSION

[t,C] = ode23s(@(t,C)Batch_Diffusion(t,C,f,kd,V0,M0,MW_M,mM0,MW_S,S0,kp0,kt0,kpD0,ktD0,Cn,Cfm,Cfs,CRD),tau,C0);
I = C(:,1);
M = C(:,2);
S = C(:,3);
X = C(:,4);
% Mass Calculations
mS = S.*MW_S.*V0;
mM = M0.*(1-X).*MW_M.*V0;
mP = mM0-mM;
wP = M0./(mP+mS+mM);

% Kinetic Constants
kp = (1/kp0+exp(Cn*wP)./kpD0).^(-1);
kt = ((1/kt0+exp(Cn*wP)./ktD0).^(-1)+CRD*kp.*(1-wP));
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
taufs = 1./(kfs.*M);

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


%% RESOLUTION WITHOUT DIFFUSION
Cfm_ND = Cfm;       % [-]
Cfs_ND = Cfs;
Ct_ND = Ct;         % [-]
kp_ND = kp0;        % [lit/mol/s]
kfm_ND = Cfm_ND.*kp_ND; % [lit/mol/s]
kfs_ND = Cfs_ND.*kp_ND;
kt_ND = kt0;        % [-]
ktc_ND = kt_ND./(1+Ct_ND);  % [lit/mol/s]
ktf_ND = ktc_ND.*Ct_ND;     % [lit/mol/s]
ktd_ND = ktD0;


[t_ND,C_ND] = ode23s(@(t,C)Batch_NoDiffusion(t,C,f,kd,kp0,kt0,Cfm_ND,Cfs_ND),tau,C0);
I_ND = C_ND(:,1);
M_ND = C(:,2);
S_ND = C(:,3);
X_ND = C(:,4);

% Characteristic Times
taup_ND = 1./(kp_ND.*M);
tautc_ND = 1./(ktc_ND.*R);
tautd_ND = 1./(ktd_ND.*R);
taufm_ND = 1./(kfm_ND.*M);
taufs_ND = 1./(kfs_ND.*M);

%Evaluate alpha
beta_ND = taup_ND./tautc_ND;
gamma_ND = taup_ND./taufm_ND+taup_ND./tautd_ND+taup_ND./taufs_ND;
alpha_ND = beta_ND + gamma_ND;

% Instantaneous Properties
mu_1_ND = 1./(alpha_ND+beta_ND./2);
mu_2_ND = (2.*gamma_ND+3.*beta_ND)./(alpha_ND.^2.*(gamma_ND+0.5*beta_ND));

% Polymer Instantaneous Properties
Mn_inst_av_ND = MW_M.*mu_1_ND;
Mw_inst_av_ND = MW_M.*mu_2_ND./mu_1_ND;


figure(1)
plot(tau,Mw_inst_av)
hold on 
plot(tau,Mw_inst_av_ND)
grid on; legend('With Diffusion','Without Diffusion','Location','best')

%% Final Result
Ratio = Mw_inst_av(end)/Mw_inst_av_ND(end);
disp(['The value of "R" equals to: ', num2str(Ratio)])

%% Function 1 - With Diffusion
function dF = Batch_Diffusion(t,C,f,kd,V0,M0,MW_M,mM0,MW_S,S0,kp0,kt0,kpD0,ktD0,Cn,Cfm,Cfs,CRD)
% Unknowns

I = C(1);
M = C(2);
S = C(3);
X = C(4);

% Preliminar Calculations
% Calculations
mM   = M0.*(1-X).*MW_M.*V0;     %[g], concentrazione*volume*peso molecolare=massa, concentrazione moltiplicata per (1-X) perchè non valutiamo le condizioni iniziali                
mS   = S0.*(1-X).*MW_S.*V0;     %[g]
mP   = mM0-mM;                  %[g], il polimero viene prodotto e all'inizio è 0 quindi la sua massa sarà la differenza tra la quantità di monomero iniziale e attuale
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

%% Function 2 -- Without Diffusion
function dF = Batch_NoDiffusion(t,C,f,kd,kp0,kt0,Cfm_ND,Cfs_ND)
% Unknowns

I = C(1);
M = C(2);
S = C(3);
X = C(4);

kp    = kp0;          %[L/mol/s] nel caso in cui non ci sia diffusione, le costanti rimangono uguali nel processo
kt    = kt0;          %[L/mol/s] nel caso in cui non ci sia diffusione, le costanti rimangono uguali nel processo
kfm   = kp.*Cfm_ND;   %[L/mol/s] ND=no diffusion
kfs   = kp.*Cfs_ND;   %[L/mol/s]

R = sqrt((2.*f.*kd.*I)./kt);

dI = -kd.*I;
dM = -(kp+kfm).*M.*R;
dS = -kfs.*S.*R;
dX = (kp+kfm).*(1-X).*R;

dF = [dI dM dS dX]';
end
