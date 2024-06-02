close all, clear all, clc
global k11 k1p kpp K_11 K_1p K_pp T Pre N y_H2O k_a
% Data
Meg = 18; %[gr/mol] - 
Mru =  72; %[gr/mol] - Molecular weight of Repeating Unit
P0_LA = 1; %[mol/lit] - Lactic Acid Concentration
W0 = 0.1; %[mol/lit] - Water Concentration
N = 60; % Number of repeating unit
T = 120 + 273.15; %[K]
Pre = 1*1.015e+5; %[Pa] - Reactor Pressure
k_a = 0.5; %[-]     factor related to water removal
y_H2O = 1; %[-]

% MM Evaluation
nru = 1:N;
MMi = Mru.*nru + Meg; % MM of the species with chain length i

% Kinetics Constants
% Direct Propagation Constant
k11 = exp(-14552./T+33.22); %P1+P1
k1p = exp(-14552./T+33.22); %P1+Pn (n>1)
kpp = exp(-14552./T+33.36); %Pn+Pm (n,m>1)
%Equilibrum Constants
K_11 = 0.33;    %[-]
K_1p = 0.275;    %[-]
K_pp = (K_1p.^2)/K_11;    %[-]

%% Resolution Section
tspan = 0:0.2:10;   %integration domain
% Initial Condition
P0 = [P0_LA zeros(1,N-1)]';
Y0 = [P0;W0];

% ODE Solver
[t,Y] = ode23s(@(t,Y)PLA(t,Y),tspan,Y0);

P= Y(:,1:N);
W = Y(:,end);

%Moments 
lambda0 = sum(P(:,1:N)'); % all chains present during the reactions
conversion = (lambda0(1) - lambda0)./lambda0(1);
lambda1 = sum([1:N]' .*P(:,1:N)');
lambda2 = sum([1:N]'.^2.*P(:,1:N)');

%Mn Mw PDI
Mn = Meg + Mru.*lambda1./lambda0;
Mw = Meg + Mru.*lambda2./lambda1;
PDI = Mw./Mn;
%% Figure Section

% Concentration of Pi species (Small oligomers --> low Ns)
N_plot = 16;
cc = jet(N_plot);
figure(1)
for i=1:N_plot
    plot(t,P(:,i),'Color',cc((i),:))
    hold on
    legendinfo1{i} = strcat(['P' num2str(i)]);
end
legend(legendinfo1)
xlabel('time [hr]'); ylabel('Concentration [mol/Lit]'); grid on


figure(2)
subplot(1,2,1)
plot(t,lambda0)
hold on
plot(t,lambda1)
plot(t,lambda2)
xlabel('time [hr]'); ylabel('Moment');
grid on
legend('\lambda_0','\lambda_1','\lambda_2','Location','best')
subplot(1,2,2)
plot(t,lambda1)
xlabel('Time [hr]'); ylabel('Moment'); grid on
legend('\lambda_1')

% Mn Mw and PDI
figure(3)
subplot(3,1,1)
plot(tspan,Mw);
xlabel('time'); ylabel('Mw'); grid on
subplot(3,1,2)
plot(tspan,Mn); xlabel('time'); ylabel('Mn'); grid on
subplot(3,1,3)
plot(tspan,PDI); xlabel('time'); ylabel('PDI'); grid on


%% Function Section
function F = PLA(t,Y)
global k11 k1p kpp K_11 K_1p K_pp T Pre N y_H2O k_a
% Species 
P = Y(1:N);
W = Y(N+1);

% Reaction Rates
r(1) = -2.*k11.*P(1).^2 + 2.*k11/K_11.*P(2).*W -2.*k1p.*P(1).*sum(P(2:N))+2.*k1p/K_1p.*W.*sum(P(3:N)); %P1+P1 --> P2
r(2) = k11.*P(1).^2-k11/K_11.*P(2).*W + ... % P1+P1 --> P2 + W
    -2.*k1p.*P(1).*P(2)+2.*k1p/K_1p.*W.*P(3) + ...% P1+P2 --> P3+W
    -2.*kpp.*P(2).*sum(P(2:N))+2.*kpp/K_pp.*W.*sum(P(4:N));
r(3) = 2.*k1p.*P(1).*(P(2)-P(3))+2.*k1p/K_1p.*W.*(P(4)-P(3))+ ...  % P1+P2 --> P3+W and P1+P3 --> P4+W
    -2.*kpp.*P(3).*sum(P(2:N))+ 2.*kpp/K_pp.*W.*sum(P(5:N)); % P3 + Pn ---> Pn3 + W

for n=3:N-2
    r(n) = 2.*k1p.*P(1).*(P(n-1)-P(n))+2.*k1p/K_1p.*W.*(P(n+1)-P(n)) + ... %P1+Pn-1 --> Pn+W and P1+Pn --> Pn+1 + W
        kpp.*(sum(P(2:n-2)'.*flip(P(2:n-2)')) - 2*P(n).*sum(P(2:N)) + ... % Pi+Pn-1 --> Pn + W
        1./K_pp.*W.*(2.*sum(P((n+2):N))-(n-3).*P(n))); % Pn + Pi ---> Pn+i + W (i: 1-->N)
end
r(N-1) = 0;
r(N) = 0;

rw = -sum(r); 

for n=1:N
    dPn(n) = r(n);
end
% Consider water in the vapor phase as an ideal gas, and using ideal gas EOS:
W_gas = Pre.*y_H2O./(8.314*T*1000);
phi = k_a.*(W-W_gas);
dw = rw - phi;
F = [dPn';dw];
end

