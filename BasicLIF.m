% Parameters of the LIF model

El = -75; %mV
V_th = -50; %mV
vReset = -80; %mV 
Rm = 100; %MΩ
Cm = 100; %pF
Ek = -80; %mV
GSRA = 1; %nS
tSRA = 200; %ms
Iapp = 500;
%V = El;
T = 150;
dt = 0.1;
t = 0:dt:T;
%GSRA = 0

%EL = -75mV, Vth = -50mV, Vreset = -80mV, Rm =100ΜΩ, Cm = 100pF, Ek = -80mV, ΔGSRA = 1 nS, and τSRA = 200ms.
%Initially set V = EL and GSRA = 0

% Initialize membrane potential
V = El * ones(1, length(t)); 
spike_train = zeros(length(t)); % To track spikes

% Simulating the LIF neuron
for i = 2:length(t)
    % Update membrane potential using Euler method
    dV = ((El - V / Rm) + GSRA(Ek - V) + Iapp)/dt;
    V(i) = V(i-1) + dV;
    
    % Check for spike
    if V(i) >= V_th
        V(i) = V_reset;   % Reset the potential
        spike_train(i) = 1; % Record a spike
    end
end

% Plotting the results
figure;
subplot(2,1,1);
plot(t, V, 'LineWidth', 2);
xlabel('Time (ms)','FontSize',14);
ylabel('Membrane Potential (mV)','FontSize',14);
title('Leaky Integrate-and-Fire Model','FontSize',24);
ylim([-80 -40]);

% Plotting spike train
subplot(2,1,2);
stem(t, spike_train, 'Marker', 'none');
xlabel('Time (ms)','FontSize',14);
ylabel('Spike','FontSize',14);
title('Spike Train','FontSize',24);
xlim([0 T]);

