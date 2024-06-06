%% DATA
% Initial Conditions
clear all, clc
beta = 6;      % You Can Change it!!!!
Min = 3; % [mol/L]
Iin = 0.0075*(1+beta); % [mol/L]
Xin = 0; % [-]
f = 0.5; % [-]

% Kinetic Rate Coefficients

kd = 2e-4; % [s-1]
Eap = 3.43*1e+4; % [J/mol]
Ap = 1.88*1e7; % [L/mol/s]
dVp = -27*1e-3; % [L/mol]
Eat = 4.6*1e3; % [J/mol]
At = 1.6*1e9; % [L/mol/s]
dVt = 15.6*1e-3; % [L/mol]

%% Resolution
global R
T = [20+beta, 100+beta]'+273.15;     %[K]
P = 1.5*1e+5; % [Pa]
R = 8.314;
% Solver
% Integration domain
tau = [1:0.1:4*3600]';        %[s]
C0 = [Iin Min Xin]';

for i=1:length(T)
    [t,C] = ode23s(@(t,C)PFR(t,C,f,kd,Eap,Ap,dVp,Eat,At,dVt,P,T(i)),tau,C0);
    Initiator(:,i) = C(:,1);
    Monomer(:,i) = C(:,2);
    Conversion(:,i) = C(:,3);
end

M_T1 = Monomer(end,1);
M_T2 = Monomer(end,2);
Ratio = M_T1/M_T2;
disp(['The R ratio is equal to: ', num2str(Ratio)])

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