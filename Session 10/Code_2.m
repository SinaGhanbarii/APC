%% FREE RADICAL POLYMERIZATION (FRP) IN CSTR AND PFR REACTORS - part 2
% Applied Physical Chemistry (APC)
% AA (2023-2024)

clear all; close all; clc
%% DATA
% Initial Conditions
% global Min Iin Xin f kd Eap Ap dVp Eat At dVt 
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

%% Resolution
global R
T = [25:10:150]' + 273.15;     %[K]
P = 101325; % [Pa]
R = 8.314;
% Solver
% Integration domain
tau = [1:0.1:1e+4]';        %[s]

% Initial Conditions
C0 = [Iin Min Xin]';        %[mol/lit]

for i=1:length(T)
    [t,C] = ode23s(@(t,C)PFR(t,C,f,kd,Eap,Ap,dVp, Eat,At,dVt,P,T(i)),tau,C0);

    I_PMA(:,i) = C(:,1);
    M_PMA(:,i) = C(:,2);
    X_PMA(:,i) = C(:,3);
end

%% Figure
cc = turbo(length(T));
figure(1)
for i=1:length(T)
    plot(t,I_PMA(:,i),'Color',cc(i,:))
    hold on
    legendinfo1{i} = strcat(['Temperature ' num2str(T(i))]);
end
xlabel('Tau [s]'); ylabel('Concentration [mol/lit]'); grid on ; title('Initiator Plot')
legend(legendinfo1)

figure(2)
for i=1:length(T)
    plot(t,M_PMA(:,i),'Color',cc(i,:))
    hold on
    legendinfo1{i} = strcat(['Temperature ' num2str(T(i))]);
end
xlabel('Tau [s]'); ylabel('Concentration [mol/lit]'); grid on ; title('Monomer Plot')
legend(legendinfo1)

figure(3)
for i=1:length(T)
    plot(t,X_PMA(:,i),'Color',cc(i,:))
    hold on
    legendinfo1{i} = strcat(['Temperature ' num2str(T(i))]);
end
xlabel('Tau [s]'); ylabel('Conversion [-]'); grid on ; title('Conversion Plot')
legend(legendinfo1)

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