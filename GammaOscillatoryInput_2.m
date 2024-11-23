% Clear previous runs and environment
clf;
clear;
close all;
clc;

% Time parameters
dt = 0.1;                    % Time step (ms)
t = 0:dt:1000;               % Time vector (ms)

% LIF Neuron parameters
Cm = 0.01;                      % Membrane capacitance (uF/cm^2)
gL = 0.002;                  % Further reduce leak conductance (mS/cm^2)
EL = -70;                    % Resting potential (mV)
Vth = -55;                   % Threshold potential (mV)
Vreset = -70;                % Reset potential (mV)
tau_m = Cm / gL;             % Membrane time constant (ms)
V = EL * ones(size(t));      % Initialize membrane potential
refractory_period = 2;        % Refractory period duration (ms)
refractory_counter = 0;       % Counter for tracking refractory period

% Input current parameters (oscillatory input)
I0 = 30;                     % Further increase amplitude of oscillatory input (nA)
f_gamma = 30;                % Frequency of gamma oscillation (Hz)
noise_amplitude = 10;          % Significantly increase noise amplitude (nA)
I_noise = noise_amplitude * randn(size(t));  % Gaussian white noise

% Generate oscillatory input (gamma frequency range) + noise
I_input = I0 * sin(2 * pi * f_gamma * t / 1000) + I_noise;  % Gamma oscillation

% Initialize spike times and firing rate tracking
spike_train = zeros(size(t));  % Spike train
V_trace = V;                   % Store voltage for plotting

% Simulation loop
for i = 2:length(t)
    if refractory_counter > 0
        % If in refractory period, hold the membrane potential at reset value
        V(i) = Vreset;
        refractory_counter = refractory_counter - dt;
    else
        % Update membrane potential using Euler method
        dV = (-(V(i-1) - EL) + I_input(i-1)) / tau_m * dt;
        V(i) = V(i-1) + dV;

        % Check if the membrane potential reaches the threshold
        if V(i) >= Vth
            V(i) = 50;                 % Set spike value for visualization
            spike_train(i) = 1;        % Mark the time of the spike
            refractory_counter = refractory_period;  % Set refractory counter
            
            % Hold the membrane potential at spike value for a few steps
            if i+10 <= length(t)       % Ensure index doesn't exceed array size
                V(i+1:i+10) = 50;      % Hold the spike for 1 ms (10 steps)
            end
            
            V(i+11) = Vreset;          % Reset membrane potential
        end
    end
    V_trace(i) = V(i);  % Record the voltage trace for plotting
end

% Plot the membrane potential and input current over time
figure;

% Plot membrane potential
subplot(3,1,1);
plot(t, V_trace, 'k', 'LineWidth', 1);
xlabel('Time (ms)');
ylabel('Membrane Potential (mV)');
title('LIF Neuron Membrane Potential with Gamma Oscillatory Input');
ylim([Vreset-5 60]);
grid on;

% Plot the input current (gamma oscillation + noise)
subplot(3,1,2);
plot(t, I_input, 'b', 'LineWidth', 1);
xlabel('Time (ms)');
ylabel('Input Current (nA)');
title('Noisy Gamma Oscillatory Input Current');
grid on;

% Plot the spike train (raster plot)
subplot(3,1,3);
plot(t, spike_train, 'k.', 'MarkerSize', 10);
xlabel('Time (ms)');
ylabel('Spikes');
title('Spike Train');
ylim([0 1.5]);
grid on;
