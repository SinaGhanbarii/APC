%% FREE RADICAL POLYMERIZATION (FRP) IN CSTR AND PFR REACTORS - part 2
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

% Kinetic Rate Coefficients

kd = 2e-4; % [s-1]
Eap = 17.7*1e3; % [J/mol]
Ap = 1.66.*1e7; % [L/mol/s]
dVp = -11.7.*1e-3; % [L/mol]
Eat = 6.7*1e3; % [J/mol]
At = 6*1e9; % [L/mol/s]
dVt = 20*1e-3; % [L/mol]