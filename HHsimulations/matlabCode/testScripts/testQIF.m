% QIF Neuron Simulation in MATLAB
% Parameters
clearvars;
close all;

Tmax = 1000;
dt = 0.01;
t = 0:dt:Tmax;

a = 0.05;
alpha = 4;
eps = 0.01;
lda = 0;
Iapp = 2;

Vth = 30;
Vrst = -65;
wrst = 4;

V = zeros(1,length(t));
w = zeros(1,length(t));
V(1) = Vrst;
w(1) = wrst;

for j=1:length(t)-1
    k1v = a*(V(j))^2-w(j)+Iapp;
    k1w = eps*(alpha*V(j)-lda-w(j));
    av = V(j)+k1v*dt;
    aw = w(j)+k1w*dt;
    k2v = av^2-aw+Iapp;
    k2w = eps*(alpha*av-lda-aw);
    V(j+1) = V(j)+(k1v+k2v)*dt/2;
    w(j+1) = w(j)+(k1w+k2w)*dt/2;
    if V(j+1) > Vth
        V(j) = 10;
        V(j+1) = Vrst;
        w(j+1) = wrst;
    end
end

figure
hold on
plot(t,V,'-b','linewidth',2);
set(gca,'fontsize',24)
% axis([0 Tmax -10 10])
xlabel('t');
ylabel('V')
