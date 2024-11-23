% Clear previous runs and environment
clf;
clear;
close all;
clc;

% Time parameters
dt = 0.1;                    % Time step (ms)
t = 0:dt:1000;               % Time vector (ms)

% LIF Neuron parameters
Cm = 1;                      % Membrane capacitance (uF/cm^2)
gL = 0.01;                   % Leak conductance (mS/cm^2)
EL = -70;                    % Resting potential (mV)
Vth = -55;                   % Threshold potential (mV)
Vreset = -70;                % Reset potential (mV)
tau_m = Cm / gL;             % Membrane time constant (ms)
V = EL * ones(size(t));      % Initialize membrane potential
refractory_period = 2;        % Refractory period duration (ms)
refractory_counter = 0;       % Counter for tracking refractory period

% Input current parameters (constant input for testing)
I_input = 200 * ones(size(t));  % Constant input current of 200 nA for testing

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

       % Update the spike handling in the simulation loop:
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

% Plot the membrane potential over time
figure;
subplot(2,1,1);
plot(t, V_trace, 'k', 'LineWidth', 1);
xlabel('Time (ms)');
ylabel('Membrane Potential (mV)');
title('LIF Neuron Membrane Potential with Constant Input');
ylim([Vreset-5 60]);
grid on;

% Plot the spike train (raster plot)
subplot(2,1,2);
plot(t, spike_train, 'k.', 'MarkerSize', 10);
xlabel('Time (ms)');
ylabel('Spikes');
title('Spike Train');
ylim([0 1.5]);
grid on;
