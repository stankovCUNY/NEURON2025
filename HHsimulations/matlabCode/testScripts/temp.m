function temp
    
    time = 0.5; % in minutes

    % Constants
    Vrest     = -65; % mV 
    dt        = 1/30; % ms
    totalTime = time*1e3; % ms
    
    C_m  = 1.0;     % membrane capacitance, in uF/cm^2

    g_Na = 120.0;   % maximum conductances, in mS/cm^2
    g_K  = 36.0;
    g_L  = 0.3;

    E_Na   = 115 + Vrest; % mV
    E_K    = -6 + Vrest; %mV
    E_L    = 10.6 + Vrest; % mV
    
    % Integration
    % The time to integrate over
    t = 0:dt:totalTime;
%     V = EIFEulerMethod(t,dt);
%     V = LIFeulerMethod(t,dt);
%     V = AdExEulerMethod(t,dt);
    
    % Initial conditions
    X0 = [Vrest; 0.053; 0.6; 0.32];

    [V,m,h,n] = HHeulerMethod(t, X0, dt);

    % Extract variables from the results
%     V = X(:, 1);
%     m = X(:, 2);
%     h = X(:, 3);
%     n = X(:, 4);

    % Calculate currents
    ina = g_Na * m.^3 .* h .* (V - E_Na);
    ik  = g_K  * n.^4 .* (V - E_K);
    il  = g_L  *         (V - E_L);

    % Plot the results
%     plot(V)
    
%     figure;
%     subplot(4,1,1);
%     
    plot(t, V, 'k');
%     ylabel('V (mV)');
%     title('Hodgkin-Huxley Neuron');
%     subplot(4,1,2);
%     plot(t, ina, 'c', t, ik, 'y', t, il, 'm');
%     ylabel('Current');
%     legend('$I_{Na}$', '$I_{K}$', '$I_{L}$');
%     subplot(4,1,3);
%     plot(t, m, 'r', t, h, 'g', t, n, 'b');
%     ylabel('Gating Value');
%     legend('m', 'h', 'n');
%     subplot(4,1,4);
%     plot(t, I_inj(t), 'k');
%     xlabel('t (ms)');
%     ylabel('$I_{inj}$ ($\mu{A}/cm^2$)');
%     ylim([-1, 31]);

end

% Define a function for the system of ODEs
function [V,m,h,n] = HHeulerMethod(t, X, dt)
    
    Vrest = -65; % mV

    C_m  = 1.0;     % membrane capacitance, in uF/cm^2

    g_Na = 120.0;   % maximum conductances, in mS/cm^2
    g_K  = 36.0;
    g_L  = 0.3;

    E_Na   = 115  + Vrest; % mV
    E_K    = -6   + Vrest; %mV
    E_L    = 10.6 + Vrest; % mV

    V = X(1);
    m = X(2);
    h = X(3);
    n = X(4);
    
    for i = 1:length(t) - 1
        
        % Calculate membrane potential & activation variables
        I_Na = g_Na * m(i)^3 * h(i) * (V(i) - E_Na);
        I_K  = g_K  * n(i)^4 *        (V(i) - E_K);
        I_L  = g_L  *                 (V(i) - E_L);

        V(i + 1) = V(i) + ((I_inj(t(i)) - I_Na - I_K - I_L) / C_m)*dt;
        m(i + 1) = m(i) + (alpha_m(V(i)) * (1.0 - m(i)) - beta_m(V(i)) * m(i))*dt;
        h(i + 1) = h(i) + (alpha_h(V(i)) * (1.0 - h(i)) - beta_h(V(i)) * h(i))*dt;
        n(i + 1) = n(i) + (alpha_n(V(i)) * (1.0 - n(i)) - beta_n(V(i)) * n(i))*dt;
        
    end
    
end


% Channel gating kinetics
% Functions of membrane voltage
% equations from original paper (Izhikevich textbook)

function alpha_n_out = alpha_n(V)
    Vrest       = -65; % mV 
    V           =  V - Vrest;
    alpha_n_out = (0.01*(10-V))/(exp((10-V)/10)-1);
end

function beta_n_out = beta_n(V)
    Vrest       = -65; % mV
    V           =  V - Vrest;
    beta_n_out  = 0.125*exp(-V/80);
end

function alpha_m_out = alpha_m(V)
    Vrest       = -65; % mV
    V           =  V - Vrest;
    alpha_m_out = (0.1*(25-V))/(exp((25-V)/10)-1);
end

function beta_m_out = beta_m(V)
    Vrest       = -65; % mV
    V           =  V - Vrest;
    beta_m_out  = 4*exp(-V/18);
end

function alpha_h_out = alpha_h(V)
    Vrest       = -65; % mV
    V           =  V - Vrest;
    alpha_h_out = 0.07*exp(-V/20);
end

function beta_h_out = beta_h(V)
    Vrest       = -65; % mV
    V           =  V - Vrest;
    beta_h_out  = 1/(exp((30-V)/10)+1);
end

%% BRETTE & GERSTNER 2005

function V = AdExEulerMethod(t,dt)
    
    % constants from table 1
    C_m     = 281;
    g_L     = 30;
    E_L     = -70.6;
    V_T     = -50.4;
    delta_T = 2;
    tau_w   = 144;
    a       = 4;
    b       = 0.0805;
    
    V(1) = E_L;
    w(1) = -1000;
    
    for i = 1:length(t) - 1
        I_L     = g_L*(V(i) - E_L);
        I_adapt = g_L*delta_T*exp((V(i)-V_T)/delta_T) - w(i);
        
        V(i + 1) = V(i) + ((I_inj(t(i)) - I_L + I_adapt) / C_m)*dt;
        w(i + 1) = w(i) + ((a*(V(i) - E_L)-w(i))/tau_w)*dt;
        
        % spike!!!
        if V(i + 1) > 30 % mV
           
            V(i + 1) = E_L;
            w(i + 1) = w(i) + b;
            
        end
        
    end
end

function V = EIFEulerMethod(t,dt)
    
    % constants from table 1
    C_m     = 281;
    g_L     = 30;
    E_L     = -70.6;
    V_T     = -50.4;
    delta_T = 2;
    
    V(1) = E_L;
    
    for i = 1:length(t) - 1
        I_L     = g_L*(V(i) - E_L);
        I_exp   = g_L*delta_T*exp((V(i)-V_T)/delta_T);
        
        V(i + 1) = V(i) + ((I_inj(t(i)) - I_L + I_exp) / C_m)*dt;
        
        % spike!!!
        if V(i + 1) > 30 % mV
            V(i + 1) = E_L;
        end
        
    end
end

function V = LIFeulerMethod(t,dt)
    
    % constants copied from HH function
    
    Vrest   = -65; % mV 
    C_m     = 1;
    g_L     = 0.3;
    E_L     = 10.6 + Vrest;
    E_K     = -6   + Vrest; %mV
    
    V(1) = Vrest;
    
    for i = 1:length(t) - 1
        I_L     = g_L*(V(i) - E_L);
        
        V(i + 1) = V(i) + ((I_inj(t(i)) - I_L) / C_m)*dt;
        
        % spike!!!
        if V(i + 1) > 30 % mV
           V(i + 1) = E_K;
        end
        
    end
end

%%

% injected current function
function I_inj_out = I_inj(t)
    if t < 50
        I_inj_out = 0;
    elseif t >= 50 && t < 60
        I_inj_out = 50; 
    elseif t >= 60
        I_inj_out = 0;
    end
end


