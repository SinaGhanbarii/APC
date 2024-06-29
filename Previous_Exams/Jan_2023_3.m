clear
clc
close all

beta = 8;
M_mass = 10; %g/L !!!!Attenzione!!!
MWm = 104.15; %g/mol 
M = M_mass/MWm; % mol/L
I = 0.005;
f = 0.8;
kd = 4e-5;
kp = 3.9e5;
kt = 6e7;
kfs = 5e3;
S = (7.8+beta)*10^-4;  % mol/ L !!!!!! Pay attention to the mothafucker concentrations' units!!!!
% if Termination only by disproportionation ---> Dpn = nu_av
% if Termination only by combination ---> Dpn = 2*nu_av
% 80% combination and 20% disproportionation ---> Dpn = 0.7*2*nu_av + 0.3*nu_av = 1.8nu_av


combination_percentage = 70;
disproportionation_percentage = 30;

X_value = 2*(combination_percentage/2)/((combination_percentage/2) + disproportionation_percentage) + disproportionation_percentage/((combination_percentage/2) + disproportionation_percentage); 


R = sqrt(2*f*kd*I/kt);

%% Without chain transfer

landa_without_chain = kp*M/(kt*R);

%% With chain transfer 

% complete formula considering monomer and initiator: landa_chain = kp*M/(kt*R + kfs*S + kfm*M + kfCTA*CTA);

landa_chain = kp*M/(kt*R+kfs*S);


%% DPn

Dpn_a = X_value*landa_without_chain;
Dpn_b = X_value*landa_chain;

R = Dpn_a/Dpn_b;
disp(['The value of R equals to: ',num2str(R)])

