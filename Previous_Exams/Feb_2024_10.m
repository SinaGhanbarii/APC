clear all; close all; clc
%% DATA
% Initial Conditions
alpha = 6;      % You Can Change it!!!!
Min = 5; % [mol/L]
Iin = 6e-3; % [mol/L]
Xin = 0; % [-]
f = 0.5; % [-]

% Kinetic Rate Coefficients

kd = 2e-4; % [s-1]
Eap = 17.71*1e3; % [J/mol]
Ap = 1.66.*1e7; % [L/mol/s]
dVp = -11.7.*1e-3; % [L/mol]
Eat = 6.7*1e3; % [J/mol]
At = 6*1e9; % [L/mol/s]
dVt = 20*1e-3; % [L/mol]

%% Resolution
global R
T = [20+10*alpha, 120+10*alpha]'+273.15;     %[K]
P = 101325; % [Pa]
R = 8.314;
% Solver
% Integration domain
tau = [1:0.1:1e+4]';        %[s]

% Initial Conditions
C0 = [Iin Min Xin]';        %[mol/lit]

for i =1:length(T)
    [t, C] = ode23s(@(t,C)PFR(t,C,f,kd,Eap,Ap,dVp,Eat,At,dVt,P,T(i)),tau,C0);
    time = t;
    Initiator(:,i) = C(:,1);
    Monomer(:,i) = C(:,2);
    Conversion(:,i) = C(:,3);
end

figure(1)
for i=1:length(T)
    subplot(1,length(T),i)
    plot(tau,Conversion(:,i))
    xlabel('Tau [s]'); ylabel('Conversion [-]'); grid on; title(['for Temperature', num2str(T(i)), 'K'])
end
% subplot(1,2,1)
% plot(tau,Conversion(:,1))
% xlabel('Tau [s]'); ylabel('Conversion [-]'); grid on; title(['for Temperature', num2str(T(1)), 'K'])
% subplot(1,2,2)
% plot(tau,Conversion(:,2))
% xlabel('Tau [s]'); ylabel('Conversion [-]'); grid on; title(['for Temperature', num2str(T(2)), 'K'])

for i=1:length(T)
    index = find(Conversion(:,i)<0.8001 & Conversion(:,i)>0.7999);
    disp(['The residence time for temperature', num2str(T(i)), 'is: ', num2str(index(end)), 's.'])
end

%% Function 
% Function 1
function F = k(T,Ea,A,dV,P)
global R
F = A.*exp(-(Ea+dV.*P*1e-6)./R./T);
end

% Function 2
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