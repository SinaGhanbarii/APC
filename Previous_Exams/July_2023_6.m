% Polymerization of styrene monomer in CSTR reactor
clear all; close all; clc

%% DATA
% Initial Conditions
beta = 6;
I0 = 0.2; % [mol/L]
f = 0.45; % [-]
T = 50 + 273.15; % [K]
P = 2e5; % [Pa]
Ap = 4.27e7; % [L/mol/s]
Ep = 3.25e4; % [J/mol]
deltaVp = -12.1e-3; % [L/mol]
At = 2.2e9; % [L/mol/s]
Et = 6.5e3; % [J/mol]
deltaVt = 14e-3; % [L/mol]

tau = 6*3600; % [s]

% Calculate kp and kt
kp = Ap*exp(-((Ep + deltaVp*P*1e-6)/(8.314*T))); % [L/mol/s]
kt = At*exp(-((Et + deltaVt*P*1e-6)/(8.314*T))); % [L/mol/s]

% Target conversion
X_target = 0.4 + 0.01*beta; % Target conversion


% Solve for kd using fzero
kd = fzero(@(kd) conversion_equation(kd, kp, kt, tau, f, I0, X_target), 1e-4);

% Display result
fprintf('The value of kd that ensures %.2f%% conversion is %.6e [s^-1]\n', X_target*100, kd);

% Verify the result
X_CSTR = 1 - 1 / (1 + kp * tau * sqrt(2*f*kd*I0/kt));
fprintf('Verification: Achieved conversion is %.2f%%\n', X_CSTR*100);

% Function to solve for kd
function F = conversion_equation(kd, kp, kt, tau, f, I0, X_target)
    X = 1 - 1 / (1 + kp * tau * sqrt(2*f*kd*I0/kt));
    F = X - X_target;
end