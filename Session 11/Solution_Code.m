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

