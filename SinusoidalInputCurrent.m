% Clear the environment and close figures
clf;
close all;
clear;
clc;

T=10;
% Define time and sinusoidal input parameters
t = 0:0.1:T;  % Time vector (1000 data points)
I0 = 5;          % Amplitude of input current (5 nA)
f = 8;           % Frequency in Hz (1 cycle per second)
phi = 0;         % Phase offset

% Calculate the sinusoidal input current
I = I0 * sin(2 * pi * f * t + phi);  % Sinusoidal function

% Create a new figure and hold it for plotting
figure;
hold on;
xlabel('Time (ms)');
ylabel('Input Current (nA)');
title('Sinusoidal Input Current to Neuron');

% Plot the sinusoidal input using plot
plot(t, I, 'b', 'LineWidth', 1.5);  % 'b' for blue color

% Set axis limits to ensure proper scaling
xlim([0 T]);  % Limit x-axis from 0 to 100 ms
ylim([-I0 I0]);  % Limit y-axis from -5 to 5 nA

% Enable grid for better visualization
grid on;
