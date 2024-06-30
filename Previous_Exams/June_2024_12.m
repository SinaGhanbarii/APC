clear, close, clc
% JUST CHANGE ALPHA!
alpha = 7;

%% DATA
Min = 5; % [mol/L]
Iin = 5e-3; % [mol/L]
Xin = 0; % [-]
f = 0.5; % [-]
kd = 2*1e-4;    %[1/s]
kp = 20000 + 1000*alpha;    %[L/mol/s]
kt = 5.1*1e+8;      %[L/mol/s]

%% RESOLUTION -- PFR
target_conversion = 0.9;
tau_min = 0;
tau_max = 10000;
tolerance = 1e-6;

while (tau_max - tau_min) > tolerance
    tau_mid = (tau_min + tau_max) / 2;
    [~, C] = ode23s(@(t,C)PFR(t,C,f,kd,kp,kt), [0 tau_mid], [Iin Min Xin]);
    X = C(end, 3);
    
    if X < target_conversion
        tau_min = tau_mid;
    else
        tau_max = tau_mid;
    end
end

tau_result = tau_mid;

%% CSTR Calculations
N_CSTR = 2+alpha;

% Single CSTR
X_CSTR_single = 1-1./((1+kp.*tau_result.*sqrt(2*f.*kd.*Iin./kt)));  % [-]

% Series of CSTRs
tau_CSTR_series = tau_result./N_CSTR;
X_CSTR_series = 1-1./((1+kp.*tau_CSTR_series.*sqrt(2*f.*kd.*Iin./kt)).^(N_CSTR));  %[-]

%% RESULTS
fprintf('Exact residence time for Conversion 0.9 in PFR: %.6f s\n', tau_result)
fprintf('Conversion in a single CSTR: %.6f\n', X_CSTR_single)
fprintf('Conversion in a series of %d CSTRs: %.6f\n', N_CSTR, X_CSTR_series)

%% FUNCTION
function dF = PFR(t,C,f,kd,kp,kt)
% Species
I = C(1);
M = C(2);
X = C(3);

% Preliminary Calculations
R = sqrt((2.*f.*kd.*I)./kt);% Overall radical concentration
% Reaction rates / mass balances
dI = -kd.*I;
dM = -kp.*R.*M;
dX = kp.*R.*(1-X);  %dX = (-1/Min*dM) --> -1/Min*(-kp*R*Min(1-X))
dF = [dI dM dX]';
end