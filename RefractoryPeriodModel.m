% Parameters of the LIF model with Refractory Period
tau_m = 20;        % Membrane time constant (ms)
R = 4;             % Membrane resistance (MΩ)
V_rest = -65;      % Resting potential (mV)
V_th = -50;        % Threshold potential (mV)
V_reset = -70;     % Reset potential after spike (mV)
I = 5;           % Input current (nA)
T = 100;           % Total simulation time (ms)
dt = 0.1;          % Time step (ms)
t = 0:dt:T;        % Time vector

% New parameter: Refractory period
refractory_period = 10;    % Refractory period (ms)
ref_time = 0;             % Keeps track of time after spike

% Initialize membrane potential
V = V_rest * ones(1, length(t)); 
spike_train = zeros(1, length(t)); % To track spikes

% Simulating the LIF neuron with refractory period
for i = 2:length(t)
    if ref_time > 0  % Check if in refractory period
        ref_time = ref_time - dt;  % Decrease the refractory time
        V(i) = V_reset;            % Hold the membrane at reset value
    else
        % Update membrane potential using Euler method
        dV = (-(V(i-1) - V_rest) + R * I) * (dt / tau_m);
        V(i) = V(i-1) + dV;
        
        % Check for spike
        if V(i) >= V_th
            V(i) = V_reset;        % Reset the potential
            spike_train(i) = 1;    % Record a spike
            ref_time = refractory_period; % Enter refractory period
        end
    end
end

% Plotting the results
figure;
subplot(2,1,1);
plot(t, V, 'LineWidth', 2);
xlabel('Time (ms)');
ylabel('Membrane Potential (mV)');
title('Leaky Integrate-and-Fire Model with Refractory Period');
ylim([-80 0]);

% Plotting spike train
subplot(2,1,2);
stem(t, spike_train, 'Marker', 'none');
xlabel('Time (ms)');
ylabel('Spike');
title('Spike Train');
xlim([0 T]);
