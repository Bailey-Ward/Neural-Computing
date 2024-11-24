% Parameters for the Hodgkin-Huxley model
Cm = 1.0;          % Membrane capacitance (uF/cm^2)
GL = 0.12;          % Leak conductance (mS/cm^2)
GNa = 120;         % Sodium channel conductance (mS/cm^2)
GK = 36;           % Potassium channel conductance (mS/cm^2)
EL = -60;      % Leak reversal potential (mV)
ENa = 45;          % Sodium reversal potential (mV)
EK = -77;          % Potassium reversal potential (mV)

% Applied current (can be a function of time)
Iapp = @(t) 10 * (t > 100 && t < 200); % Example: current applied between t=10ms and t=50ms

% Time parameters
tspan = [0 100];   % Simulation time in ms
dt = 0.01;         % Time step in ms
time = tspan(1):dt:tspan(2);

% Initial conditions
Vm = EL; % Initial membrane potential (mV)
m = 0; % Initial sodium activation gating variable
h = 0;  % Initial sodium inactivation gating variable
n = 0; % Initial potassium activation gating variable

% Preallocate arrays for results
Vm_array = zeros(size(time));
m_array = zeros(size(time));
h_array = zeros(size(time));
n_array = zeros(size(time));

Vm_array(1) = Vm;
m_array(1) = m;
h_array(1) = h;
n_array(1) = n;

% Helper functions for gating variable dynamics
%alpha_m = @(V) 0.1 * (V + 40) ./ (1 - exp(-(V + 40) / 10));
%beta_m = @(V) 4 * exp(-(V + 65) / 18);

%alpha_h = @(V) 0.07 * exp(-(V + 65) / 20);
%beta_h = @(V) 1 ./ (1 + exp(-(V + 35) / 10));

%alpha_n = @(V) 0.01 * (V + 55) ./ (1 - exp(-(V + 55) / 10));
%beta_n = @(V) 0.125 * exp(-(V + 65) / 80);

alpha_m =@(Vm) 10^5 * (-Vm - 0.045) / (exp(100 * (-Vm- 0.045)) - 1);
beta_m =@(Vm) 4 * 10^3 * exp((-Vm - 0.070) / 0.018);

alpha_h =@(Vm) 70 * exp(50 * (-Vm - 0.070));
beta_h =@(Vm) 10^3 / (1 + exp(100 * (-Vm - 0.040)));

alpha_n =@(Vm) (10^4 * (-Vm - 0.060)) / (exp(100 * (-Vm - 0.060)) - 1);
beta_n =@(Vm) 125 * exp((-Vm - 0.070) / 0.08);


% Euler method for numerical integration
for i = 1:length(time)-1
    % Update gating variables using forward Euler
    m = m + dt * (alpha_m(Vm) * (1 - m) - beta_m(Vm) * m);
    h = h + dt * (alpha_h(Vm) * (1 - h) - beta_h(Vm) * h);
    n = n + dt * (alpha_n(Vm) * (1 - n) - beta_n(Vm) * n);
    
    % Compute ionic currents
    INa = GNa * m^3 * h * (ENa - Vm);
    IK = GK * n^4 * (EK - Vm);
    IL = GL * (EL - Vm);
    
    % Update membrane potential
    dVm = (INa + IK + IL + Iapp(time(i))) / Cm;
    Vm = Vm + dt * dVm;
    
    % Store results
    Vm_array(i+1) = Vm;
    m_array(i+1) = m;
    h_array(i+1) = h;
    n_array(i+1) = n;
end

% Plot results
figure;
subplot(2,1,1);
plot(time, Vm_array, 'LineWidth', 1.5);
xlabel('Time (ms)');
ylabel('Membrane Potential (mV)');
title('Membrane Potential Dynamics');
grid on;

subplot(2,1,2);
plot(time, m_array, 'b', 'LineWidth', 1.5); hold on;
plot(time, h_array, 'r', 'LineWidth', 1.5);
plot(time, n_array, 'g', 'LineWidth', 1.5);
xlabel('Time (ms)');
ylabel('Gating Variables');
legend('m', 'h', 'n');
title('Gating Variable Dynamics');
grid on;
