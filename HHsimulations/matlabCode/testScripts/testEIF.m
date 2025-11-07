% Time parameters
dt = 1/30;               % Time step (ms)
T = 100;                 % Total simulation time (ms)
time = 0:dt:T;           % Time vector

% Neuron parameters
C = 281;                 % Membrane capacitance (pF)
gL = 30;                 % Leak conductance (nS)
EL = -60;              % Leak reversal potential (mV)
DeltaT = 2;              % Slope factor (mV)
VT = -50.4;              % Spike threshold (mV)
V_reset = -65;         % Reset potential (mV)

% Input current
Iext = 300;              % External current (pA)

% Initialization
V = EL;                  % Membrane potential (mV)
V_trace = zeros(size(time));  % Voltage trace

for t = 1:length(time)

    dVdt = (-gL*(V - EL) + gL*DeltaT*exp((V - VT)/DeltaT) + Iext)/C;
    V = V + dVdt * dt;

    if V >= VT
%         V = 60;
%     elseif V == 60
        V = V_reset;
    end

    V_trace(t) = V;
end

% Plot the results
figure;
plot(time, V_trace);
xlabel('Time (ms)');
ylabel('Membrane Potential (mV)');
title('Exponential Integrate-and-Fire Neuron Model');