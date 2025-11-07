tic

% Parameters of the Izhikevich model
a = 0.02;
b = -1;
c = -55;
d = 6;

% Simulation time parameters
T = 1*60*1000; % Total time in ms
dt = 1/30; % Time step in ms
n = T/dt; % Number of iterations

% External input current
I = 50; % Constant input current

% Initial conditions
v = -65; % Initial membrane potential
u = 0; % Initial recovery variable

% Preallocate vectors for speed
v_values = zeros(1, n);
u_values = zeros(1, n);

for t = 1:n
    % Update variables using Izhikevich model equations
    dv = 0.04 * v^2 + 5 * v + 140 - u + I;
    du = a * (b * v - u);
    
    v = v + dv * dt;
    u = u + du * dt;
    
    % Check for spike and apply reset conditions
    if v >= 30
        v = c;
        u = u + d;
    end
    
    % Store values for plotting
    v_values(t) = v;
    u_values(t) = u;
end

% Plot the results
figure;
subplot(3,1,1);
plot((0:1e4-1) * dt, v_values(1:1e4));
title('Membrane potential (v) over time');
xlabel('Time (ms)');
ylabel('v (mV)');

subplot(3,1,2);
plot((0:1e4-1) * dt, u_values(1:1e4));
title('Recovery variable (u) over time');
xlabel('Time (ms)');
ylabel('u');

subplot(3,1,3);
histogram(v_values)
title('v histogram');
xlabel('v (mV))');
ylabel('Counts');

toc