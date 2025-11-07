% plot(wave4Ana)
% hold on
% gain = (1/300e3)/(150e-12);
% waveIntegral = -cumtrapz((1:72)*(1/30000),wave4Ana)*gain;
% plot(waveIntegral)
% % plot(diff(waveIntegral))
% hold off

%%

g_cleft   = 1e-6;
g_ephapse = 1e-22;
C_ephapse = 1e-12;

% C1 = 1;
% C2 = 1;

C1 = g_ephapse/C_ephapse;
C2 = (g_cleft - g_ephapse)/C_ephapse;

gt = 1:72;
g  = wave4Ana-wave4Ana(1);
g  = g*1e-3; % convert to volts

tspan = [1 72];
ic = 0;
opts = odeset('RelTol',1e-2,'AbsTol',1e-4);
[t,y] = ode45(@(t,y) myode(t,y,gt,g,C1,C2), tspan, ic, opts);

figure
plot(1:72,g)
hold on

F = griddedInterpolant(t,y);

plot(1:72,F(1:72)-g')
hold off

%%

vt = wave4Ana;
vt = vt(28:51);
vt = vt - vt(1);

tiledlayout(2,1)

nexttile
plot((-8:15)*(1/30),vt)

nexttile
plot((-8:15)*(1/30),cumtrapz((g_ephapse/C_ephapse)*vt))
 

