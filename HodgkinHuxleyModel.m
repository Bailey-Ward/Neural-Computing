% Hodgkin-Huxley Neuron Model Simulation
clear all

% Parameters
C_m = 100;    % Membrane capacitance, uF/cm^2
g_Na = 12.0; % Maximum conductance of sodium (mS/cm^2)
g_K = 3.6;   % Maximum conductance of potassium (mS/cm^2)
g_L = 30;    % Leak conductance (mS/cm^2)
E_Na = 45;    % Sodium reversal potential (mV)
E_K = -82;    % Potassium reversal potential (mV)
E_L = -60;  % Leak reversal potential (mV)
V_rest = -65; % Resting membrane potential (mV)

% Time parameters
dt = 0.01;  % Time step (ms)
T = 350;    % Total time (ms)
time = 0:dt:T;  % Time vector

% External current (stimulation)
I_ext = zeros(size(time));
I_ext(500:600) = 10;  % Inject current between 5 and 6 ms (in microamps)

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
alpha_m = @(V) (0.1 * (V + 40)) / (1 - exp(-(V + 40) / 10));
beta_m = @(V) 4.0 * exp(-(V + 65) / 18);
alpha_h = @(V) 0.07 * exp(-(V + 65) / 20);
beta_h = @(V) 1.0 / (1 + exp(-(V + 35) / 10));
alpha_n = @(V) (0.01 * (V + 55)) / (1 - exp(-(V + 55) / 10));
beta_n = @(V) 0.125 * exp(-(V + 65) / 80);

% Simulation loop
for t = 2:length(time)
    % Update gating variables using Euler method
    m = m + dt * (alpha_m(V(t-1)) * (1 - m) - beta_m(V(t-1)) * m);
    h = h + dt * (alpha_h(V(t-1)) * (1 - h) - beta_h(V(t-1)) * h);
    n = n + dt * (alpha_n(V(t-1)) * (1 - n) - beta_n(V(t-1)) * n);
    
    % Compute conductances for sodium and potassium
    g_Na_t = g_Na * (m^3) * h;
    g_K_t = g_K * (n^4);
    
    % Compute ionic currents
    I_Na = g_Na_t * (V(t-1) - E_Na);
    I_K = g_K_t * (V(t-1) - E_K);
    I_L = g_L * (V(t-1) - E_L);
    
    % Update membrane potential using Euler's method
    V(t) = V(t-1) + dt * (I_ext(t) - (I_Na + I_K + I_L)) / C_m;
    
    % Store values for plotting
    m_values(t) = m;
    h_values(t) = h;
    n_values(t) = n;
    V_values(t) = V(t);
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
