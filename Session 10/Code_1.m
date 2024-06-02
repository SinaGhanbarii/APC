%% FREE RADICAL POLYMERIZATION (FRP) IN CSTR AND PFR REACTORS - part 1
% Applied Physical Chemistry (APC)
% AA (2023-2024)

clear all
close all
clc

%% DATA
% Initial Conditions

Min = 4; % [mol/L]
Iin = 6e-3; % [mol/L]
Xin = 0; % [-]
f = 0.5; % [-]
T = 50 + 273.15; % [K]
P = 101325; % [Pa]

% Kinetic rate constants for the monomers:
% Ethylene Styrene MethylMethacrylate MethAcrylate

kd = [2 2 2 2]'.*1e-4; % [s-1]
kp = [54 238 648 22900]';
kt = [2.9 2.0 0.94 5.1].'*1e8;

%% Resolution 
% PFR and single CSTR for 4 monomers 
% Integration domain
tau = 1:0.1:1e+4;    %[s] --- tspan

% Initial Conditions
C0 = [Iin Min Xin]';

for i=1:length(kd)
    [t,C] = ode23s(@(t,C)PFR(t,C,f,kd(i),kp(i),kt(i)),tau,C0);
    %PFR Numerical Solution
    I_PFR_Numerical(:,i) = C(:,1);      % [mol/lit]
    M_PFR_Numerical(:,i) = C(:,2);      % [mol/lit]
    X_PFR_Numerical(:,i) = C(:,3);      % [mol/lit]

    %PFR Analytical Solution
    I_PFR_Analytical(:,i) = Iin.*exp(-kd(i).*tau);      % [mol/lit]
    R_PFR(:,i) = sqrt((2.*f.*kd(i).*I_PFR_Analytical(:,i))/kt(i));      % [mol/lit]
    M_PFR_Analytical(:,i) = Min.*exp(2.*kp(i)./kd(i).*sqrt(2.*f.*kd(i).*Iin./kt(i)).*(exp(-kd(i).*tau./2)-1));      % [mol/lit]
    X_PFR_Analytical(:,i) = 1-exp(2.*kp(i)./kd(i).*sqrt(2.*f.*kd(i).*Iin./kt(i)).*(exp(-kd(i).*tau./2)-1));      % [mol/lit]

    % CSTR Analytical Solution - no dF
    I_CSTR(:,i) = Iin./(1+kd(i).*tau);      % [mol/lit]
    R_CSTR(:,i) = sqrt(2.*f.*kd(i).*I_CSTR(:,i)./kt(i));      % [mol/lit]
    M_CSTR(:,i) = Min./(1+kp(i).*R_CSTR(:,i).*tau');      % [mol/lit]
    X_CSTR(:,i) = 1-(1./(1+kp(i).*R_CSTR(:,i).*tau'));      % [mol/lit]
end
%% Figures
figure(1)
% M_PFR Analytical and numerical 
subplot(1,2,1)
plot(t,M_PFR_Analytical)
hold on 
set(gca,'ColorOrderIndex')
plot(t,M_PFR_Numerical,'--','LineWidth',2)
xlabel('Tau [s]'); ylabel('Concentration of Monomer [Mol/Lit]'); grid on
legend('PE','PS','PMMA','PMA','PE_{Num}','PS_{Num}','PMMA_{Num}','PMA_{Num}')

% X_PFR Analytical and numerical 
subplot(1,2,2)
plot(t,X_PFR_Analytical)
hold on 
set(gca,'ColorOrderIndex')
plot(t,X_PFR_Numerical,'--','LineWidth',2)
xlabel('Tau [s]'); ylabel('Conversion [-]'); grid on
legend('PE','PS','PMMA','PMA','PE_{Num}','PS_{Num}','PMMA_{Num}','PMA_{Num}')


figure(2) %Comparison I_CSTR and I_PFR_Numerical
plot(t,I_PFR_Numerical)
hold on
set(gca,'ColorOrderIndex',1)
plot(t,I_CSTR,'--','LineWidth',1)
xlabel('Tau [s]'); ylabel('Concentration [mol/lit]'); grid on
legend('PE_{PFR}','PS_{PFR}','PMMA_{PFR}','PMA_{PFR}','PE_{CSTR}','PS_{CSTR}','PMMA_{CSTR}','PMA_{CSTR}')

figure(3) % Comparison M_CSTR and M_PFR
plot(t,M_PFR_Numerical)
hold on
set(gca,'ColorOrderIndex',1)
plot(t,M_CSTR,'--','LineWidth',1)
xlabel('Tau [s]'); ylabel('Concentration [mol/lit]'); grid on
legend('PE_{PFR}','PS_{PFR}','PMMA_{PFR}','PMA_{PFR}','PE_{CSTR}','PS_{CSTR}','PMMA_{CSTR}','PMA_{CSTR}')

figure(4)
plot(t,X_PFR_Numerical)
hold on
set(gca,'ColorOrderIndex',1)
plot(t,X_CSTR,'--','LineWidth',1)
xlabel('Time [s]'); ylabel('Conversion [-]'); grid on
legend('PE_{PFR}','PS_{PFR}','PMMA_{PFR}','PMA_{PFR}','PE_{CSTR}','PS_{CSTR}','PMMA_{CSTR}','PMA_{CSTR}')

%% Resolution for PMA in PFR and Series of CSTRs
N_CSTR = 10;

for j=1:N_CSTR
    tau_CSTR_series(:,j) = tau./j; %[s]

    I_CSTR_Series(:,j) = Iin./((1+kd(4).*tau_CSTR_series(:,j)).^(j));   %[mol/lit]
    M_CSTR_Series(:,j) = Min./((1+kp(4).*tau_CSTR_series(:,j).*sqrt(2*f.*kd(4).*I_CSTR_Series(:,j)./kt(4))).^(j));  %[mol/lit]
    X_CSTR_Series(:,j) = 1-1./((1+kp(4).*tau_CSTR_series(:,j).*sqrt(2*f.*kd(4).*I_CSTR_Series(:,j)./kt(4))).^(j));  %[-]
    
end

%% Figures for PMA in PFR and Series of CSTRs
cc = turbo(N_CSTR);
figure(5)   % Comparison I
for j=1:N_CSTR
    plot(tau,I_CSTR_Series(:,j),'Color',cc(j,:))
    hold on
    legendinfo1{j} = strcat(['Reactor Number ' num2str(j)]);
end
xlabel('Tau [s]'); ylabel('Concentration [mol/lit]');title('Initiator') ;grid on
legend(legendinfo1,'Location','best')

figure(6)   % Comparison M
for j=1:N_CSTR
    plot(tau,M_CSTR_Series(:,j),'Color',cc(j,:))
    hold on
    legendinfo1{j} = strcat(['Reactor Number ' num2str(j)]);
end
xlabel('Tau [s]'); ylabel('Concentration [mol/lit]');title('Monomer') ;grid on
legend(legendinfo1,'Location','best')

figure(7)   % Comparison X
for j=1:N_CSTR
    plot(tau,X_CSTR_Series(:,j),'Color',cc(j,:))
    hold on
    legendinfo1{j} = strcat(['Reactor Number ' num2str(j)]);
end
xlabel('Tau [s]'); ylabel('Conversion [-]');title('Conversion') ;grid on
legend(legendinfo1,'Location','best')
%% Function

function dF = PFR(t,C,f,kd,kp,kt)
% Species
I = C(1);
M = C(2);
X = C(3);

% Preliminary Calculations
R = sqrt((2.*f.*kd.*I)./kt);% Overall radical concentration
% Reaction rates / mass balances
dI = -kd.*I;
dM = -kp.*R.*M;
dX = kp.*R.*(1-X);  %dX = (-1/Min*dM) --> -1/Min*(-kp*R*Min(1-X))
dF = [dI dM dX]';
end