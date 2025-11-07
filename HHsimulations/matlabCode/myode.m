function dydt = myode(t,y,gt,g,C1,C2)
g = interp1(gt,g,t); % Interpolate the data set (gt,g) at time t
dydt = -C1*y + C2*g; % Evaluate ODE at time t
