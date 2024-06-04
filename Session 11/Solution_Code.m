%% Session 11 APC
clear all, close all, clc

%% Data  
V0 = 20;        % [L]
M0 = 4.7;       % [mol/lit]
w_I0 = 0.01;    % [w/w]
MW_M = 100;     % [gr/mol]
MW_S = 78;      % [gr/mol]
MW_I = 164;     % [gr/mol]
rho_M = 940;    % [gr/lit]
rho_S = 880;    % [gr/lit]
f = 0.5;        % [-]
kd = 5.55e-6;   % [s^-1]
kp0 = 715;      % [lit/mol/s]
kt0 = 9.8e6;    % [lit/mol/s]
kpD0 = 3e+11;   % [lit/mol/s]
ktD0 = 3e8;     % [lit/mol/s]
Cn = 25;        % [-]
CRD = 180;      % [-]

%% Resolution
% Preliminar Calculations
mM0 = M0*V0*MW_M;   % [gr]
mS0 = (V0-mM0/rho_M)*rho_S; % [gr]
mI0 = w_I0*mM0;     % [gr]
I0 = (mI0/MW_I)/V0; % [mol/lit]

% Target Conversions
X_target = [0.1 0.4 0.7 0.9]';
% Initial Conditions
C0 = [I0 M0]';      % [mol/lit]
n = 1:1:1e+4;
%% Diffusion - point a) b) c)
Cfm = 0;    Ct = 1e+3; % Case 1A
% Cfm = 0;    Ct = 0.001; % Case 1B
% Cfm = 0.01; Ct = 1e+3;  % Case 1C

Cfm_ND = Cfm;       % [-]
Ct_ND = Ct;         % [-]
kp_ND = kp0;        % [lit/mol/s]
kfm_ND = Cfm_ND.*kp_ND; % [lit/mol/s]
kt_ND = kt0;        % [-]
ktc_ND = kt_ND./(1+Ct_ND);  % [lit/mol/s]
ktf_ND = ktc_ND.*Ct_ND;     % [lit/mol/s]
ktd_ND = ktD0;
% Solver
for i=1:length(X_target)
    % Integration Domain Definition
    X_span(:,i) = 0:(X_target(i)/100):X_target(i);

    % Diffusion
    [X,C] = ode23s(@(X,C)Batch_Diffusion(X,C,f,kd,V0,M0,MW_M,mM0,mS0,kp0,kt0,kpD0,ktD0,Cn,Cfm,CRD),X_span(:,i),C0);
    I(:,i) = C(:,1);        % [mol/lit]
    M(:,i) = C(:,2);        % [mol/lit]
    
    % Calculations
    mM(:,i) = M0.*(1-X_target(i)).*MW_M.*V0;    % [gr]
    mP(:,i) = mM0-mM(:,i);  % [gr]
    wp(:,i) = mP(:,i)./(mP(:,i)+mS0+mM(:,i));       % [gr/gr]

    % Kinetic Constants
    kp(:,i) = (1./kp0+exp(Cn.*wp(:,i))./kpD0).^(-1);    % [lit/mol/s]
    kt(:,i) = ((1./kt0+exp(Cn.*wp(:,i))./ktD0).^(-1))+CRD*kp(:,i).*(1-wp(:,i));    % [lit/mol/s]
    kfm(:,i) = kp(:,i).*Cfm;            % [lit/mol/s]
    ktc(:,i) = kt(:,i)./(1+Ct);         % [lit/mol/s]
    ktd(:,i) = Ct.*kt(:,i)./(1+Ct);     % [lit/mol/s]  ---> kt = ktc + ktd

    % Characteristic times
    tau_p(:,i) = 1./(kp(:,i).*M(end,i));    % Propagation Characteristic Time
    tau_fm(:,i) = 1./(kfm(:,i).*M(end,i));  % CTM Characteristic Time 

    R(:,i) = sqrt(2.*f.*kd.*I(end,i)./kt(:,i));     % [mol/lit]
    
    tau_tc(:,i) = 1./(ktc(:,i).*R(:,i));    % Termination By Combination Characteristic Time
    tau_td(:,i) = 1./(ktd(:,i).*R(:,i));    % Termination By Propogation Characteristic Time

    % alpha, beta and gamma
    beta(:,i) = tau_p(:,i)./tau_tc(:,i);    % [-]
    gamma(:,i) = tau_p(:,i)./tau_fm(:,i)+tau_p(:,i)./tau_td(:,i);   %[-]
    alpha(:,i) = beta(:,i)+gamma(:,i);      % [-]

    % Moments 
    mu_1(:,i) = 1./(gamma(:,i)+beta(:,i)./2);       % alpha << 1
    % mu_2(:,i) = (2.*(gamma(:,i)+3.*beta(:,i))./(alpha(:,i).^(2).*gamma(:,i)+0.5.*beta(:,i)));
    mu_2(:,i) = (2*gamma(:,i)+3*beta(:,i))./(alpha(:,i).^(2).*(gamma(:,i)+0.5*beta(:,i)));

    % Instantaneous number and length CLD
    xn_inst(:,i) = (alpha(:,i)./(1+alpha(:,i)).^(n)).*((gamma(:,i)...
        +0.5.*(n-1).*alpha(:,i).*beta(:,i))./(gamma(:,i)+0.5.*beta(:,i)));  % [-]
    xw_inst(:,i) = n.*(alpha(:,i)./(1+alpha(:,i)).^(n)).*((gamma(:,i)...
        +0.5.*(n-1).*alpha(:,i).*beta(:,i))./(gamma(:,i)+0.5.*beta(:,i)))./mu_1(:,i);   % [-]
    % Instantaneous MWs, PDI, and DPs
    Mn_av_inst(:,i) = mu_1(:,i).*MW_M;  % [gr/mol]
    Mw_av_inst(:,i) = mu_2(:,i)./mu_1(:,i).*MW_M;   % [gr/mol]
    PDI_inst(:,i) = Mw_av_inst(:,i)./Mn_av_inst(:,i);
    
    DPn_inst(:,i) = Mn_av_inst(:,i)./MW_M;  % [-]
    DPw_inst(:,i) = Mw_av_inst(:,i)./MW_M;  % [-]


    %% No Diffusion
    [X_ND, C_ND] = ode23s(@(X,C)Batch_NoDiffusion(X,C,f,kd,kp0,kt0,Cfm_ND),X_span(:,i),C0);
    I_ND(:,i) = C_ND(:,1);        % [mol/lit]
    M_ND(:,i) = C_ND(:,2);        % [mol/lit]
    
    % Characteristic times
    tau_p_ND(:,i) = 1./(kp_ND.*M_ND(end,i));    % Propagation Characteristic Time
    tau_fm_ND(:,i) = 1./(kfm_ND.*M_ND(end,i));  % CTM Characteristic Time 

    R_ND(:,i) = sqrt(2.*f.*kd.*I_ND(end,i)./kt_ND);     % [mol/lit]
    
    tau_tc_ND(:,i) = 1./(ktc_ND.*R_ND(:,i));    % Termination By Combination Characteristic Time
    tau_td_ND(:,i) = 1./(ktd_ND.*R_ND(:,i));    % Termination By Propogation Characteristic Time

    % alpha, beta and gamma
    beta_ND(:,i) = tau_p_ND(:,i)./tau_tc_ND(:,i);    % [-]
    gamma_ND(:,i) = tau_p_ND(:,i)./tau_fm_ND(:,i)+tau_p_ND(:,i)./tau_td_ND(:,i);   %[-]
    alpha_ND(:,i) = beta_ND(:,i)+gamma_ND(:,i);      % [-]

    % Moments 
    mu_1_ND(:,i) = 1./(gamma_ND(:,i)+beta_ND(:,i)./2);       % alpha << 1
    % mu_2_ND(:,i) = (2.*(gamma_ND(:,i)+3.*beta_ND(:,i))./(alpha_ND(:,i).^(2).*gamma_ND(:,i)+0.5.*beta_ND(:,i)));
    mu_2_ND(:,i) = (2.*gamma_ND(:,i)+3.*beta_ND(:,i))./(alpha_ND(:,i).^(2).*(gamma_ND(:,i)+0.5*beta_ND(:,i)));


    % Instantaneous number and length CLD
    xn_inst_ND(:,i) = (alpha_ND(:,i)./(1+alpha_ND(:,i)).^(n)).*((gamma_ND(:,i)...
        +0.5.*(n-1).*alpha_ND(:,i).*beta_ND(:,i))./(gamma_ND(:,i)+0.5.*beta_ND(:,i)));  % [-]
    xw_inst_ND(:,i) = n.*(alpha_ND(:,i)./(1+alpha_ND(:,i)).^(n)).*((gamma_ND(:,i)...
        +0.5.*(n-1).*alpha_ND(:,i).*beta_ND(:,i))./(gamma_ND(:,i)+0.5.*beta_ND(:,i)))./mu_1_ND(:,i);   % [-]
    % Instantaneous MWs, PDI, and DPs
    Mn_av_inst_ND(:,i) = mu_1_ND(:,i).*MW_M;  % [gr/mol]
    Mw_av_inst_ND(:,i) = mu_2_ND(:,i)./mu_1_ND(:,i).*MW_M;   % [gr/mol]
    PDI_inst_ND(:,i) = Mw_av_inst_ND(:,i)./Mn_av_inst_ND(:,i);
    
    DPn_inst_ND(:,i) = Mn_av_inst_ND(:,i)./MW_M;  % [-]
    DPw_inst_ND(:,i) = Mw_av_inst_ND(:,i)./MW_M;  % [-]

end
%% Figures --- All Cases
figure(1)
subplot(1,4,1)
hold on 
plot(n,xn_inst,'linewidth',2); 
xlabel('Chain length n [-]')
ylabel('xn Inst [-]')
title('xn Inst')
leg1 = legend( 'X =0.1', 'X =0.4', 'X =0.7', 'X =0.9');

subplot(1,4,2)
hold on 
plot(n,xw_inst,'linewidth',2); 
xlabel('Chain length n [-]')
ylabel('xw Inst [-]')
title('xw Inst')
leg1 = legend( 'X =0.1', 'X =0.4', 'X =0.7', 'X =0.9');

subplot(1,4,3)
hold on 
plot(n,xn_inst_ND,'linewidth',2); 
xlabel('Chain length n [-]')
ylabel('xn Inst ND [-]')
title('xn Inst ND')
leg1 = legend( 'X =0.1', 'X =0.4', 'X =0.7', 'X =0.9');

subplot(1,4,4)
hold on 
plot(n,xw_inst_ND,'linewidth',2); 
xlabel('Chain length n [-]')
ylabel('xw Inst ND [-]')
title('xw Inst ND')
leg1 = legend( 'X =0.1', 'X =0.4', 'X =0.7', 'X =0.9');

figure(2) 

subplot(1,4,1)
hold on 
plot(X_target, DPn_inst, 'b--o','Linestyle','--','linewidth',2); % b (blue) / o (marker), it could be *
xlabel('Conversion [-]') 
ylabel('DPn Inst [-]') 
title('DPn Inst') 
leg1 = legend('DPn Inst');

subplot(1,4,2)
hold on 
plot(X_target, DPw_inst,'r--o','Linestyle','--','linewidth',2); 
xlabel('Conversion [-]') 
ylabel('DPw Inst [-]') 
title('DPw Inst') 
leg1 = legend('DPw Inst');

subplot(1,4,3)
hold on 
plot(X_target, DPn_inst_ND, 'b--o','Linestyle','--','linewidth',2); % b (blue) / o (marker), it could be *
xlabel('Conversion [-]') 
ylabel('DPn Inst ND[-]') 
title('DPn Inst ND') 
leg1 = legend('DPn Inst ND');

subplot(1,4,4)
hold on 
plot(X_target, DPw_inst_ND,'r--o','Linestyle','--','linewidth',2); 
xlabel('Conversion [-]') 
ylabel('DPw Inst ND [-]') 
title('DPw Inst ND')
leg1 = legend('DPw Inst ND');

figure(3)
subplot(1,4,1)
hold on 
plot(X_target, Mn_av_inst,'b--o','Linestyle','--','linewidth',2);
xlabel('Conversion [-]') 
ylabel('Mn Av Inst [g/mol]') 
title('Mn Av Inst')
leg1 = legend('Mn Inst');

subplot(1,4,2)
hold on 
plot(X_target, Mw_av_inst,'r--o','Linestyle','--','linewidth',2);
xlabel('Conversion [-]') 
ylabel('Mw Av Inst [g/mol]') 
title('Mw Av Inst')
leg1 = legend('Mw Inst');

subplot(1,4,3)
hold on 
plot(X_target, Mn_av_inst_ND,'b--o','Linestyle','--','linewidth',2);
xlabel('Conversion [-]') 
ylabel('Mn Av Inst ND [g/mol]') 
title('Mn Av Inst ND')
leg1 = legend('Mn Inst ND');

subplot(1,4,4)
hold on 
plot(X_target, Mw_av_inst_ND,'r--o','Linestyle','--','linewidth',2);
xlabel('Conversion [-]') 
ylabel('Mw Av Inst ND [g/mol]') 
title('Mw Av Inst ND')
leg1 = legend('Mw Inst ND');


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

dI = -(kd.*I)./((kp+kfm).*(1-X).*sqrt(2.*f.*kd.*I./kt));
dM = -M./(1-X);

dF = [dI dM]';
end
%% Function 2 --- No Diffusion Case
function dF = Batch_NoDiffusion(X,C,f,kd,kp0,kt0,Cfm_ND)
% Unknowns 
I = C(1);
M = C(2);

% Kinetics Constants
kp = kp0;           % [lit/mol/s]
kfm = kp.*Cfm_ND;   % [lit/mol/s]
kt = kt0;           % [lit/mol/s]

% Mass Balance
dI = -(kd.*I)./((kp+kfm).*(1-X).*sqrt(2.*f.*kd.*I./kt));
dM = -M./(1-X);

dF = [dI dM]';
end
