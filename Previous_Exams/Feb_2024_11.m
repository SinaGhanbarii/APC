%% Excercise - 11
clear all; close all; clc
alpha = 6;
% DATA
V0 = 100;
M0 = 10;
w0_I = 0.02;
MW_M = 100;
MW_S = 78;
MW_I = 164;
rho_M = 940;    % [gr/lit]
rho_S = 880;    % [gr/lit]
f = 0.5;        
kd = 5.55e-6;
kp0 = 715;
kt0 = 9.8e+6;
kpD0 = 3e+11;
ktD0 = 3e+8;
Cn = 25;
CRD = 180;

%% Resolution
% Preliminar Calculations
mM0 = M0*V0*MW_M;   % [gr]
mS0 = (V0-mM0/rho_M)*rho_S; % [gr]
mI0 = w0_I*mM0;     % [gr]
I0 = (mI0/MW_I)/V0; % [mol/lit]
X_target = 0.7+0.01*alpha;
X_span = [0:0.001:X_target]';
%%  Case 1 
Cfm_c1 = 0.0001; Ct_c1 = 900;   % Case 1
C0 = [I0 M0]';
[X_1, C_1] = ode23s(@(X,C)Batch_Diffusion(X,C,f,kd,V0,M0,MW_M,mM0,mS0,kp0,kt0,kpD0,ktD0,Cn,Cfm_c1,Ct_c1,CRD),X_span,C0);
Conversion_1 = X_1;
I_1 = C_1(:,1);
M_1 = C_1(:,2);

% Calculations
mM_1 = M0.*(1-X_target).*MW_M.*V0;    % [gr]
mP_1 = mM0-mM_1;  % [gr]
wp_1 = mP_1./(mP_1+mS0+mM_1);       % [gr/gr]

% Kinetic Constants
kp_1 = (1./kp0+exp(Cn.*wp_1)./kpD0).^(-1);    % [lit/mol/s]
kt_1 = ((1./kt0+exp(Cn.*wp_1)./ktD0).^(-1))+CRD*kp_1.*(1-wp_1);    % [lit/mol/s]
kfm_1 = kp_1.*Cfm_c1;            % [lit/mol/s]
ktc_1 = kt_1./(1+Ct_c1);         % [lit/mol/s]
ktd_1 = Ct_c1.*kt_1./(1+Ct_c1);     % [lit/mol/s]  ---> kt = ktc + ktd

% Characteristic times
tau_p_1 = 1./(kp_1.*M_1);    % Propagation Characteristic Time
tau_fm_1 = 1./(kfm_1.*M_1);  % CTM Characteristic Time 

R_1 = sqrt(2.*f.*kd.*I_1./kt_1);     % [mol/lit]
    
tau_tc_1 = 1./(ktc_1.*R_1);    % Termination By Combination Characteristic Time
tau_td_1 = 1./(ktd_1.*R_1);    % Termination By Propogation Characteristic Time

% alpha, beta and gamma
beta_1 = tau_p_1./tau_tc_1;    % [-]
gamma_1 = tau_p_1./tau_fm_1+tau_p_1./tau_td_1;   %[-]
alpha_1 = beta_1+gamma_1;      % [-]

% Moments 
mu_1_c1 = 1./(gamma_1+beta_1./2);       % alpha << 1
mu_2_c1 = (2*gamma_1+3*beta_1)./(alpha_1.^(2).*(gamma_1+0.5*beta_1));

Mn_av_inst_1 = mu_1_c1.*MW_M;  % [gr/mol]
Mw_av_inst_1 = mu_2_c1./mu_1_c1.*MW_M;   % [gr/mol]    
DPn_inst_1 = Mn_av_inst_1./MW_M;  % [-]
DPw_inst_1 = Mw_av_inst_1./MW_M;  % [-]

result_case1 = [Mn_av_inst_1(end) Mw_av_inst_1(end) DPn_inst_1(end) DPw_inst_1(end)];
table(result_case1)

%% Case 2
Cfm_c2 = 0;    Ct_c2 = 0.0009;  % Case 2
[X_2, C_2] = ode23s(@(X,C)Batch_Diffusion(X,C,f,kd,V0,M0,MW_M,mM0,mS0,kp0,kt0,kpD0,ktD0,Cn,Cfm_c2,Ct_c2,CRD),X_span,C0);
I_2 = C_2(:,1);
M_2 = C_2(:,2);

% Calculations
mM_2 = M0.*(1-X_target).*MW_M.*V0;    % [gr]
mP_2 = mM0-mM_2;  % [gr]
wp_2 = mP_2./(mP_2+mS0+mM_2);       % [gr/gr]

% Kinetic Constants
kp_2 = (1./kp0+exp(Cn.*wp_2)./kpD0).^(-1);    % [lit/mol/s]
kt_2 = ((1./kt0+exp(Cn.*wp_2)./ktD0).^(-1))+CRD*kp_2.*(1-wp_2);    % [lit/mol/s]
kfm_2 = kp_2.*Cfm_c2;            % [lit/mol/s]
ktc_2 = kt_2./(1+Ct_c2);         % [lit/mol/s]
ktd_2 = Ct_c1.*kt_2./(1+Ct_c2);     % [lit/mol/s]  ---> kt = ktc + ktd

% Characteristic times
tau_p_2 = 1./(kp_2.*M_2);    % Propagation Characteristic Time
tau_fm_2 = 1./(kfm_2.*M_2);  % CTM Characteristic Time 

R_2 = sqrt(2.*f.*kd.*I_2./kt_2);     % [mol/lit]
    
tau_tc_2 = 1./(ktc_2.*R_2);    % Termination By Combination Characteristic Time
tau_td_2 = 1./(ktd_2.*R_2);    % Termination By Propogation Characteristic Time

% alpha, beta and gamma
beta_2 = tau_p_2./tau_tc_2;    % [-]
gamma_2 = tau_p_2./tau_fm_2+tau_p_2./tau_td_2;   %[-]
alpha_2 = beta_2+gamma_2;      % [-]

% Moments 
mu_1_c2 = 1./(gamma_2+beta_2./2);       % alpha << 1
mu_2_c2 = (2*gamma_2+3*beta_2)./(alpha_2.^(2).*(gamma_2+0.5*beta_2));

Mn_av_inst_2 = mu_1_c2.*MW_M;  % [gr/mol]
Mw_av_inst_2 = mu_2_c2./mu_1_c2.*MW_M;   % [gr/mol]    
DPn_inst_2 = Mn_av_inst_2./MW_M;  % [-]
DPw_inst_2 = Mw_av_inst_2./MW_M;  % [-]

result_case2 = [Mn_av_inst_2(end) Mw_av_inst_2(end) DPn_inst_2(end) DPw_inst_2(end)];
table(result_case2)
%% Function
function dF = Batch_Diffusion(X,C,f,kd,V0,M0,MW_M,mM0,mS0,kp0,kt0,kpD0,ktD0,Cn,Cfm,Ct,CRD)
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

dI = -(kd.*I)./((kp+kfm).*(1-X).*sqrt(2.*f.*kd.*I./kt));
dM = -M./(1-X);

dF = [dI dM]';
end