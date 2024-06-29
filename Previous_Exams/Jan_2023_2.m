% In Case of Series CSTR
clear all
close all
clc
%% DATA
% Initial Conditions
beta = 8;
I0 = 0.001; % [mol/L]
f = 0.6; % [-]
T = 60 + 273.15; % [K]
P = 3e5; % [Pa]
n = 4;
Ap = 1.79e7; % L/mol/s
Ep = (17+beta)*10^3; %J
deltaVp = -11.7e-3;

kd1 = 2e-4; % [s-1]
kd2 = 2*kd1;
kp = Ap*exp(-((Ep + deltaVp*P*10^-6)/8.314/T)); % [L/mol/s]
kt = 5e15; % [L/mol/s]

tau =0.5; % [s]
j=4;
tau_CSTR = tau./j;
X_CSTR_a = 1-1./((1+kp.*tau_CSTR.*sqrt(2*f.*kd1.*I0./kt)).^(j));
X_CSTR_b = 1-1./((1+kp.*tau_CSTR.*sqrt(2*f.*kd2.*I0./kt)).^(j));


R = X_CSTR_b/X_CSTR_a;
disp(['The value of R equals to: ', num2str(R)])



