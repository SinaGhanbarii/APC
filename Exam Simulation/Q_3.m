% Excercise 3 Bulk FRP
clear all, close all, clc
%% DATA
Min = 5;
Iin = 4.8e-3;
Xin = 0;
f = 0.6;

% Kinetics Constants
kd = 2e-4;
kp = 757;
kt = 1.9e+7;

%% Resolution 
C0 = [Iin Min Xin]';
% tau = [0:1e-6:1.7e-3]';
tau = [0:0.1:1e+4]';
[t,C] = ode23s(@(t,C)Batch_Reactor(t,C,f,kd,kp,kt),tau,C0);

figure(1)
plot(tau,C(:,3)); xlabel('Tau [s]'); ylabel('Conversion [-]'); grid on
Conversion = C(:,3);
index = find(Conversion > 0.5999 & Conversion < 0.6000);
Conversion_desired = Conversion(index(end))
desired_time = tau(index(end))
% % Analytical Solution
% I_Batch_Analytical = Iin.*exp(-kd.*tau);
% R_SS = sqrt(2*f*kd.*Iin./kt).*exp(-kd.*tau./2);
% M_Batch_Analytical = Min.*exp(2.*kp./kd*sqrt(2*f*kd.*Iin./kt).*(exp(-kd.*tau./2)-1));
% X_Batch_Analytical = 1 - exp(2.*kp./kd.*sqrt(2*f*kd*Iin./kt).*(exp(-kd.*tau./2)-1));
% 
% figure(2)
% plot(tau,X_Batch_Analytical); xlabel('Tau [s]'); ylabel('Conversion [-]'); grid on


%% Figure

%% Function
function dF = Batch_Reactor(t,C,f,kd,kp,kt)
% Unknowns 
I = C(1);
M = C(2);
X = C(3);

% Mass Balances
R = sqrt((2*f*kd.*I)./kt);

dI = -kd.*I;
dM = - kp.*M.*R;
dX = kp.*R.*(1-X);

dF = [dI dM dX]';

end