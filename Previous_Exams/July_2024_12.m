%% DATA
% Initial Conditions
clear all, clc
alpha = 7;      % You can change this if needed
Min = 5; % [mol/L]
Iin = 6e-3; % [mol/L]
Xin = 0; % [-]
f = 0.5; % [-]

% Kinetic Rate Coefficients
kd = 2e-4; % [s-1]
Eap = 17.71e+3; % [J/mol]
Ap = 1.66*1e7; % [L/mol/s]
dVp = -11.7e-3; % [L/mol]
Eat = 6.7*1e3; % [J/mol]
At = 6*1e9; % [L/mol/s]
dVt = 20*1e-3; % [L/mol]

%% Resolution
global R
T = [20+10*alpha, 120+10*alpha]'+273.15;     %[K] Modified to match the question
P = 1.01325e5; % [Pa] Atmospheric pressure
R = 8.314;

% Integration settings
tspan = [0:0.1:10*3600]; % Increased upper limit to ensure X=0.8 is reached
C0 = [Iin Min Xin]';
options = odeset('Events', @(t,C) eventFunction(t,C));

residence_times = zeros(1,2);

for i=1:length(T)
    [t,C,te,Ce,ie] = ode15s(@(t,C)PFR(t,C,f,kd,Eap,Ap,dVp,Eat,At,dVt,P,T(i)),tspan,C0,options);
    residence_times(i) = te;
end

Ratio = residence_times(2)/residence_times(1);
disp(['The ratio of residence times (T2/T1) is: ', num2str(Ratio)])

%% Functions
function F = k(T,Ea,A,dV,P)
global R
F = A.*exp(-(Ea+dV.*P*1e-6)./R./T);
end

function dF = PFR(t,C,f,kd,Eap,Ap,dVp,Eat,At,dVt,P,T)
I = C(1);
M = C(2);
X = C(3);

kp = k(T,Eap,Ap,dVp,P);
kt = k(T,Eat,At,dVt,P);

R = sqrt((2.*f.*kd.*I)./kt);

% mass balances 
dI = -kd.*I;
dM = -kp.*M.*R;
dX = kp.*R.*(1-X);

dF = [dI dM dX]';
end

function [value,isterminal,direction] = eventFunction(t,C)
value = C(3) - 0.8;  % We want to stop when X reaches 0.8
isterminal = 1;  % Stop the integration
direction = 0;   % Approach from either direction
end