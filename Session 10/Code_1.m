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
end
%% Figures

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