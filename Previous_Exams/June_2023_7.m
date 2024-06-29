clear, clc
global R Ap Ep dVp T P f kd Iin Min Xin beta target_conversion

%% Input Data
Min = 1;  % Initial monomer concentration (assumed)
Iin = 0.01;  % Initial Initiator Concentration
Xin = 0;  % Initial conversion
f = 0.5;  
R = 8.314;  % [J/mol/K]
beta = 9;  % Change it!!!!!!

% Kinetics Constants
kd = 2e-4;  

Ap = 2.5e+6;  
Ep = 2.1e+4;  
dVp = -16e-3;  

T = 50+273.15;  % [K] 
P = 1e+5;  % [Pa] 

% Target conversion after 1 hour
target_conversion = 10*(beta+0.5);
target_time = 3600;  % 1 hour in seconds

% Use fzero to find kt that gives the target conversion
kt = fzero(@calculate_kt, 1e+7);

% Simulate with the calculated kt
[t, C] = ode45(@(t,C) Batch_Reactor(t, C, kt), 0:60:3600, [Iin Min Xin]');

% Display results
fprintf('Beta value: %.2f\n', beta);
fprintf('Target conversion: %.2f%%\n', target_conversion);
fprintf('Calculated kt: %.4e\n', kt);
fprintf('Final conversion: %.2f%%\n', C(end,3)*100);

% Helper functions
function F = k(T, Ea, A, dV, P)
    global R
    F = A*exp(-(Ea+dV*P*1e-6)/(R*T));
end

function dF = Batch_Reactor(t, C, kt)
    global Ap Ep dVp T P f kd
    
    I = C(1);
    M = C(2);
    X = C(3);
    
    kp = k(T, Ep, Ap, dVp, P);
    R = sqrt((2*f*kd*I)/kt);
    
    dI = -kd*I;
    dM = -kp*M*R;
    dX = kp*R*(1-X);
    
    dF = [dI; dM; dX];
end

function diff = calculate_kt(kt)
    global Iin Min Xin target_conversion
    
    [~, C] = ode45(@(t,C) Batch_Reactor(t, C, kt), [0 3600], [Iin Min Xin]);
    final_conversion = C(end, 3) * 100;  % Convert to percentage
    diff = final_conversion - target_conversion;
end