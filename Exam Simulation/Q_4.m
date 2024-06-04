clear all; close all; clc
% DATA
V0 = 50;
M0 = 8;
w0_I = 0.005;
MW_M = 100;
MW_S = 78;
MW_I = 164;
rho_M = 940;
rho_S = 880;
f = 0.5;
kd = 5.55e-6;
kp0 = 715;
kt0 = 9.8e+6;
kpD0 = 3e+11;
ktD0 = 3e+8;
Cn = 25;
CRD = 180;

%% Resolution
% Preliminary Calculations
mM0 = M0*V0*MW_M;
mS0 = (V0-mM0/rho_M)*rho_S;
mI0 = w0_I*mM0;
I0 = (mI0/MW_I)/V0;
X0 = 0;
tau = 0:0.1:1e+4;

n = 1:15*1e+3;
Cfm = 0.01;
Ct = 1000;

% Solver
tspan = [0:1:3600]';
C0 = [I0 M0 X0]';
[t C] = ode23s(@(t,C)Reaction_Time(t,C,f,kd,V0,M0,MW_M,mM0,mS0,kp0,kt0,kpD0,ktD0,Ct,Cn,Cfm,CRD),tspan,C0);

I = C(:,1);
M = C(:,2);
X = C(:,3);

% Preliminar Calculations
mS = mS0;
mM = M0.*(1-X).*MW_M.*V0;
mP = mM0-mM;
wP = M0./(mP+mS+mM);

% Kinetic Constants
kp = (1/kp0+exp(Cn*wP)./kpD0).^(-1);
kt = ((1/kt0+exp(Cn*wP)./ktD0).^(-1)+CRD*kp.*(1-wP));
kfm = kp*Cfm;
ktc = kt./(1+Ct);
ktd = Ct.*kt./(1+Ct);

R = sqrt((2.*f.*kd.*I)./kt);

% Characteristic Times
taup = 1./(kp.*M);
tautc = 1./(ktc.*R);
tautd = 1./(ktd.*R);
taufm = 1./(kfm.*M);

beta = taup./tautc;
gamma = taup./tautd + taup./taufm;
alpha = beta+gamma;

% Instantaneous Properties
mu_1 = 1./(alpha+beta./2);
mu_2 = (2.*gamma+3.*beta)./(alpha.^2.*(gamma+0.5*beta));

xn_inst = (alpha./(1+alpha).^(n)).*((gamma+0.5*(n-1).*alpha.*beta)./(gamma+0.5*beta));
xw_inst = n.*(alpha./(1+alpha).^(n)).*((gamma+0.5*(n-1).*alpha.*beta)./(gamma+0.5*beta))./mu_1;

% Polymer Instantaneous Properties
Mn_inst_av = MW_M.*mu_1;
Mw_inst_av = MW_M.*mu_2./mu_1;
DPn_inst = Mn_inst_av./MW_M;
DPw_inst = Mw_inst_av./MW_M;
PDI_inst = Mw_inst_av./Mn_inst_av;
%% Function
function dF = Reaction_Time(t,C,f,kd,V0,M0,MW_M,mM0,mS0,kp0,kt0,kpD0,ktD0,Ct,Cn,Cfm,CRD)
% Unknowns

I = C(1);
M = C(2);
X = C(3);

% Preliminar Calculations
mS = mS0;
mM = M0.*(1-X).*MW_M.*V0;
mP = mM0-mM;
wP = M0./(mP+mS+mM);

% Kinetic Constants
kp = (1/kp0+exp(Cn*wP)./kpD0).^(-1);
kt = ((1/kt0+exp(Cn*wP)./ktD0).^(-1)+CRD*kp.*(1-wP));
kfm = kp*Cfm;
ktc = kt./(1+Ct);
ktd = Ct.*kt./(1+Ct);

R = sqrt((2.*f.*kd.*I)./kt);

% Characteristic Times
taup = 1./(kp.*M);
tautc = 1./(ktc.*R);
tautd = 1./(ktd.*R);
taufm = 1./(kfm.*M);

beta = taup./tautc;
gamma = taup./tautd + taup./taufm;
alpha = beta+gamma;

dI = -kd.*I;
dM = -(kp+kfm).*M.*R;
dX = (kp+kfm).*(1-X).*R;

% Additional equation for the polymer 
% Mass Balance on P
% dP = (kp+kfm).*M.*R.*(gamma+0.5*beta)

dF = [dI dM dX]';


end