clear all
close all
clc

%% DATA
global R
beta = 6;
n_CSTR = 4;
Min = 4; % [mol/L]
Iin = 0.001; % [mol/lit]
Xin = 0;
f = 0.6;
kd = 1e-4;  
Ap = 1.79e+7;
Eap = (17+beta)*1e+3;
dVp = -11.7e-3;
Eat = 6.7*1e3; % [J/mol]
At = 6*1e9; % [L/mol/s]
dVt = 20*1e-3; % [L/mol]


%% Resolution
T = 60 + 273.15;
P = 3e+5;
R = 8.314;
tspan = [1:0.1:1e+4]';
C0 = [Iin Min Xin]';
[t,C] = ode23s(@(t,C)PFR(t,C,f,kd,Eap,Ap,dVp,Eat,At,dVt,P,T),tspan,C0);
I_PMA = C(:,1);
M_PMA = C(:,2);
X_PMA = C(:,3);





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
