% Hodgkin-Huxley Neuron Model Simulation
clear all

% Parameters
g_L = 30;    % Leak conductance (mS/cm^2)
g_Na = 12.0; % Maximum conductance of sodium (mS/cm^2)
g_K = 3.6;   % Maximum conductance of potassium (mS/cm^2)
E_Na = 45;    % Sodium reversal potential (mV)
E_K = -82;    % Potassium reversal potential (mV)
E_L = -60;  % Leak reversal potential (mV)
C_m = 100;    % Membrane capacitance, uF/cm^2
Iapp = 0; %variable applied current
V_rest = -65; % Resting membrane potential (mV)
Vm = E_L;

% Time parameters
dt = 0.01;  % Time step (ms)
T = 350;    % Total time (ms)
time = 0:dt:T;  % Time vector
Iapp(time >= 100 & time < 200) = 0.22;

% Initialize variables
V = V_rest * ones(size(time));  % Membrane potential (mV)
m = 0.05;   % Sodium activation gating variable
h = 0.6;    % Sodium inactivation gating variable
n = 0.32;   % Potassium activation gating variable

% Storage for gating variables over time
m_values = zeros(size(time));
h_values = zeros(size(time));
n_values = zeros(size(time));
V_values = zeros(size(time));

% Define rate functions for gating variables

alpha_m =@(V) (10^5 * (-Vm - 0.045)) / (exp(100 * (-Vm- 0.045)) - 1);
beta_m =@(V) 4 * 10^3 * exp((-Vm - 0.070) / 0.018);
alpha_h =@(V) 70 * exp(50 * (-Vm - 0.070));
beta_h =@(V) 10^3 / (1 + exp(100 * (-Vm - 0.040)));
alpha_n =@(V) (10^4 * (-Vm - 0.060)) / (exp(100 * (-Vm - 0.060)) - 1);
beta_n =@(V) 125 * exp((-Vm - 0.070) / 0.08);


% Simulation loop
for t = 2:length(time)

    % Update gating variables using Euler method
    m = m + dt * (alpha_m(Vm) * (1 - m) - beta_m(Vm) * m);
    h = h + dt * (alpha_h(Vm) * (1 - h) - beta_h(Vm) * h);
    n = n + dt * (alpha_n(Vm) * (1 - n) - beta_n(Vm) * n);
    
    %dm = (alpha_m * (1-m)) - (beta_m * m);
    %dh = (alpha_h * (1-h)) - (beta_h * h);
    %dn = (alpha_n * (1-n)) - (beta_n * n);

    % Compute conductances for sodium and potassium
    g_Na_t = g_Na * (m^3) * h;
    g_K_t = g_K * (n^4);
    
    % Compute ionic currents
    I_Na = g_Na_t * (V(t-1) - E_Na);
    I_K = g_K_t * (V(t-1) - E_K);
    I_L = g_L * (V(t-1) - E_L);
    
    % Update membrane potential using Euler's method
    dVm = (I_Na + I_K + I_L + Iapp(time(t))) / Cm;
    Vm = Vm + dt * dVm;

    % Store values for plotting
    m_values(t+1) = m;
    h_values(t+1) = h;
    n_values(t+1) = n;
    V_values(t+1) = Vm;
end

% Plot membrane potential over time
figure;
subplot(2,1,1);
plot(time, V_values, 'LineWidth', 2);
title('Hodgkin-Huxley Model: Membrane Potential');
xlabel('Time (ms)');
ylabel('Membrane Potential (mV)');
grid on;

% Plot gating variables over time
subplot(2,1,2);
plot(time, m_values, 'r', 'LineWidth', 1.5); hold on;
plot(time, h_values, 'g', 'LineWidth', 1.5);
plot(time, n_values, 'b', 'LineWidth', 1.5);
legend('m (Na+ activation)', 'h (Na+ inactivation)', 'n (K+ activation)');
title('Gating Variables Over Time');
xlabel('Time (ms)');
ylabel('Gating Variable');
grid on;
