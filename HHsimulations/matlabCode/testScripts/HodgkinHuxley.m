function HodgkinHuxley
    
    close all
    
    time = 10; % in minutes
    nWaveTrials = 1000;

    noiseSc = 0; % normal noise scale
    waveSc  = 500; 

    Vrest = -65; % mV − change this to −65 if desired
    dt = 1/30; % ms
    totalTime = time*1e3; % ms
    C = 1; % uF/cm^2
    
    % constants; values based on Table 1
    E_Na = 115 + Vrest; % mV
    E_K = -6 + Vrest; %mV
    E_Leak = 10.6 + Vrest; % mV
    
    g_Na = 120; % mS/cm^2
    g_K = 36; % mS/cm^2
    g_Leak = 0.3; % mS/cm^2
    
    g_Cleft = 1;
    
    % Vector of timesteps
    t = [0:dt:totalTime];
    
    % Current input −− change this to see how different inputs affect the neuron
    I_current = ones(1,length(t))*0.0;
    I_current(50/dt:end) = 0; % Input of 3 microA/cm2 beginning at 50 ms and steady until end of time period.
    
    waveform = [0;
                0.0162680427000000;
                0.00351580409999994;
                0.0172214811000001;
                0.0849553341000000;
                0.137652669000000;
                0.124483301100000;
                0.0897622527000001;
                0.0543658521000000;
                -0.0207174219000001;
                -0.125893595400000;
                -0.218893566000000;
                -0.295704947100000;
                -0.351302323800000;
                -0.351620136600000;
                -0.283052025000000;
                -0.157575558900000;
                0.0101302830000002;
                0.179226555900000;
                0.275027251800000;
                0.258798935700000;
                0.197123389200000;
                0.180239584200000;
                0.190171234200000];
    waveform = waveSc*waveform;
    
    % 
    I_wave = ones(1,length(t))*0.0;
    
    UL = length(t) - length(waveform) - 1;
    LL = 30;
    numElements = UL-LL+1;
    r = randperm(numElements, nWaveTrials);
    waveIdx = r + LL;
    
    I_wave(waveIdx) = 1;
    
    I_wave = conv(I_wave,waveform);
    I_wave = I_wave(1:length(t));
    
    I_current = I_current - I_wave;
    
    % initializing values
    V(1) = Vrest; % membrane potential is starting at its resting state
    
    % separate functions to get the alpha and beta values
    [alphaM, betaM] = m_equations(V(1), Vrest);
    [alphaN, betaN] = n_equations(V(1), Vrest);
    [alphaH, betaH] = h_equations(V(1), Vrest);
    
    % initializing gating variables to the asymptotic values when membrane potential
    % is set to the membrane resting value based on equation 13
    m(1) = (alphaM / (alphaM + betaM));
    n(1) = (alphaN / (alphaN + betaN));
    h(1) = (alphaH / (alphaH + betaH));
    
    % repeat for time determined in totalTime , by each dt
    for i = 1:length(t)
        
        % calculate new alpha and beta based on last known membrane potenatial
        [alphaN, betaN] = n_equations(V(i), Vrest);
        [alphaM, betaM] = m_equations(V(i), Vrest);
        [alphaH, betaH] = h_equations(V(i), Vrest);
        
        % conductance variables − computed separately to show how this
        % changes with membrane potential in one of the graphs
        conductance_K(i) = g_K*(n(i)^4);
        conductance_Na(i)= g_Na*(m(i)^3)*h(i);
        
        % retrieving ionic currents
        I_Na(i) = conductance_Na(i)*(V(i)-E_Na);
        I_K(i) = conductance_K(i)*(V(i)-E_K);
        I_Leak(i) = g_Leak*(V(i)-E_Leak);
        
        % Calculating the input
        Input = I_current(i) - (I_Na(i) + I_K(i) + I_Leak(i));
        
        % Calculating the new membrane potential
        V(i+1) = V(i)*(1 + noiseSc*normrnd(0,1)) + Input* dt*(1/C);
        
        % getting new values for the gating variables
        m(i+1) = m(i) + (alphaM *(1-m(i)) - betaM * m(i))*dt;
        n(i+1) = n(i) + (alphaN *(1-n(i)) - betaN * n(i))*dt;
        h(i+1) = h(i) + (alphaH *(1-h(i)) - betaH * h(i))*dt;
    end
    
    % spike times
    thresh=0;    % or whatever level you wish
    [pks,spikeTimes] = findpeaks(V,'MinPeakHeight',thresh);
    
    tiledlayout(5,1)
    
    nexttile
    plot(t,I_current)
    
    nexttile
    plot(t,V(:,1:end-1))
    
    binSize        = 1/30;
    trialLength    = 60;
    preTrialLength = 30;
    
    spikeTimes     = spikeTimes*(1/30);
    trialTimes     = [(waveIdx - preTrialLength)*(1/30);(waveIdx + trialLength)*(1/30)];
    
    color = 'k';
    
    nexttile
    [histCounts] = PSTH(spikeTimes',trialTimes',trialLength*binSize,preTrialLength*binSize,binSize,color);
    plot((-29:1:60)/30,histCounts)
    
    nexttile
    plot((0:23)/30,waveform)
    xlim([-1 2])
    
    nexttile 
    plot(t,I_Na+I_K+I_Leak+I_current)
end

% calculate alpha m and beta m based on Table 2
function [alpha_m, beta_m] = m_equations(V, Vrest)
    alpha_m = (2.5-0.1*(V-Vrest))/(exp(2.5-0.1*(V-Vrest))-1);
    beta_m = 4*exp((Vrest-V)/18);
end

% calculate alpha n and beta n based on Table 2
function [alpha_n, beta_n] = n_equations(V, Vrest)
    alpha_n = (0.1-0.01*(V-Vrest))/(exp(1-0.1*(V-Vrest))-1);
    beta_n = 0.125*exp((Vrest-V)/80);
end

% calculate alpha h and beta h based on Table 2
function [alpha_h, beta_h] = h_equations(V, Vrest)
    alpha_h = 0.07*exp((Vrest-V)/20);
    beta_h = 1/(1+exp(3-0.1*(V-Vrest)));
end