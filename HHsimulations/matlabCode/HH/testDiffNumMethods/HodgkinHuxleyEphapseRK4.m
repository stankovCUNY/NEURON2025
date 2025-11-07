function HodgkinHuxleyEphapseRK4(timeMin,inputType)
    
    tic

%     close all
    
    time = timeMin*60;% in seconds
    
    noiseSc = 0.02; % normal noise scale

    Vrest = -65; % mV − change this to −65 if desired
    dt = 1/30; % ms
    fs  = (1/dt)*1000; % hz
    totalTime = time*1e3; % ms
    C = 1; % uF/cm^2
    
    % constants; values based on Table 1
    E_Na = 115 + Vrest; % mV
    E_K = -6 + Vrest; %mV
    E_Leak = 10.6 + Vrest; % mV
    
    g_Na = 120; % mS/cm^2
    g_K = 36; % mS/cm^2
    g_Leak = 0.3; % mS/cm^2
    
    gain = 6;
    
    % Vector of timesteps
    t = [0:dt:totalTime];
    
    % target ion channels switch
    ionSwitch = ones(1,length(t));
    ionSwitch(1:60*fs) = 0;
    
    % junction switch
    juncSwitch = ones(1,length(t));
    juncSwitch(60*fs+1:120*fs) = 0;
    
    frozenNoiseRef = noiseSc*normrnd(0,1,1,length(t));
    frozenNoiseTar = noiseSc*normrnd(0,1,1,length(t));
    
    %%
    
    % Current input −− change this to see how different inputs affect the neuron
    I_current_ref = ones(1,length(t))*0.0;
    I_current_tar = ones(1,length(t))*0.0;
    
    I_current_ref(50/dt:end) = 0; % Input of 3 microA/cm2 beginning at 50 ms and steady until end of time period.
    I_current_tar(30/dt:end) = 0; % Input of 3 microA/cm2 beginning at 50 ms and steady until end of time period.
    
    % initializing values
    Vref(1)  = Vrest; % membrane potential is starting at its resting state
    
    % separate functions to get the alpha and beta values
    [alphaMref, betaMref] = m_equations(Vref(1), Vrest);
    [alphaNref, betaNref] = n_equations(Vref(1), Vrest);
    [alphaHref, betaHref] = h_equations(Vref(1), Vrest);
    
    % initializing gating variables to the asymptotic values when membrane potential
    % is set to the membrane resting value based on equation 13
    
    mRef(1) = (alphaMref / (alphaMref + betaMref));
    nRef(1) = (alphaNref / (alphaNref + betaNref));
    hRef(1) = (alphaHref / (alphaHref + betaHref));
    
    % solve for common input spike times
    % Runge-Kutta-4 (RK4) integration loop for the reference neuron
    for i = 1:length(t)-1
        % K1
        [alphaMref, betaMref] = m_equations(Vref(i), Vrest);
        [alphaNref, betaNref] = n_equations(Vref(i), Vrest);
        [alphaHref, betaHref] = h_equations(Vref(i), Vrest);
        conductance_K_ref = g_K*(nRef(i)^4);
        conductance_Na_ref = g_Na*(mRef(i)^3)*hRef(i);
        I_Na_ref(i) = conductance_Na_ref*(Vref(i)-E_Na);
        I_K_ref(i) = conductance_K_ref*(Vref(i)-E_K);
        I_Leak_ref(i) = g_Leak*(Vref(i)-E_Leak);
        Input_ref = I_current_ref(i) - (I_Na_ref(i) + I_K_ref(i) + I_Leak_ref(i));
        k1V = (Input_ref / C);
        k1m = (alphaMref *(1-mRef(i)) - betaMref * mRef(i));
        k1n = (alphaNref *(1-nRef(i)) - betaNref * nRef(i));
        k1h = (alphaHref *(1-hRef(i)) - betaHref * hRef(i));

        % K2
        Vref_k2 = Vref(i) + 0.5 * dt * k1V;
        mRef_k2 = mRef(i) + 0.5 * dt * k1m;
        nRef_k2 = nRef(i) + 0.5 * dt * k1n;
        hRef_k2 = hRef(i) + 0.5 * dt * k1h;
        [alphaMref_k2, betaMref_k2] = m_equations(Vref_k2, Vrest);
        [alphaNref_k2, betaNref_k2] = n_equations(Vref_k2, Vrest);
        [alphaHref_k2, betaHref_k2] = h_equations(Vref_k2, Vrest);
        conductance_K_ref_k2 = g_K*(nRef_k2^4);
        conductance_Na_ref_k2 = g_Na*(mRef_k2^3)*hRef_k2;
        I_Na_ref_k2 = conductance_Na_ref_k2*(Vref_k2-E_Na);
        I_K_ref_k2 = conductance_K_ref_k2*(Vref_k2-E_K);
        I_Leak_ref_k2 = g_Leak*(Vref_k2-E_Leak);
        Input_ref_k2 = I_current_ref(i) - (I_Na_ref_k2 + I_K_ref_k2 + I_Leak_ref_k2);
        k2V = (Input_ref_k2 / C);
        k2m = (alphaMref_k2 *(1-mRef_k2) - betaMref_k2 * mRef_k2);
        k2n = (alphaNref_k2 *(1-nRef_k2) - betaNref_k2 * nRef_k2);
        k2h = (alphaHref_k2 *(1-hRef_k2) - betaHref_k2 * hRef_k2);

        % K3
        Vref_k3 = Vref(i) + 0.5 * dt * k2V;
        mRef_k3 = mRef(i) + 0.5 * dt * k2m;
        nRef_k3 = nRef(i) + 0.5 * dt * k2n;
        hRef_k3 = hRef(i) + 0.5 * dt * k2h;
        [alphaMref_k3, betaMref_k3] = m_equations(Vref_k3, Vrest);
        [alphaNref_k3, betaNref_k3] = n_equations(Vref_k3, Vrest);
        [alphaHref_k3, betaHref_k3] = h_equations(Vref_k3, Vrest);
        conductance_K_ref_k3 = g_K*(nRef_k3^4);
        conductance_Na_ref_k3 = g_Na*(mRef_k3^3)*hRef_k3;
        I_Na_ref_k3 = conductance_Na_ref_k3*(Vref_k3-E_Na);
        I_K_ref_k3 = conductance_K_ref_k3*(Vref_k3-E_K);
        I_Leak_ref_k3 = g_Leak*(Vref_k3-E_Leak);
        Input_ref_k3 = I_current_ref(i) - (I_Na_ref_k3 + I_K_ref_k3 + I_Leak_ref_k3);
        k3V = (Input_ref_k3 / C);
        k3m = (alphaMref_k3 *(1-mRef_k3) - betaMref_k3 * mRef_k3);
        k3n = (alphaNref_k3 *(1-nRef_k3) - betaNref_k3 * nRef_k3);
        k3h = (alphaHref_k3 *(1-hRef_k3) - betaHref_k3 * hRef_k3);

        % K4
        Vref_k4 = Vref(i) + dt * k3V;
        mRef_k4 = mRef(i) + dt * k3m;
        nRef_k4 = nRef(i) + dt * k3n;
        hRef_k4 = hRef(i) + dt * k3h;
        [alphaMref_k4, betaMref_k4] = m_equations(Vref_k4, Vrest);
        [alphaNref_k4, betaNref_k4] = n_equations(Vref_k4, Vrest);
        [alphaHref_k4, betaHref_k4] = h_equations(Vref_k4, Vrest);
        conductance_K_ref_k4 = g_K*(nRef_k4^4);
        conductance_Na_ref_k4 = g_Na*(mRef_k4^3)*hRef_k4;
        I_Na_ref_k4 = conductance_Na_ref_k4*(Vref_k4-E_Na);
        I_K_ref_k4 = conductance_K_ref_k4*(Vref_k4-E_K);
        I_Leak_ref_k4 = g_Leak*(Vref_k4-E_Leak);
        Input_ref_k4 = I_current_ref(i) - (I_Na_ref_k4 + I_K_ref_k4 + I_Leak_ref_k4);
        k4V = (Input_ref_k4 / C);
        k4m = (alphaMref_k4 *(1-mRef_k4) - betaMref_k4 * mRef_k4);
        k4n = (alphaNref_k4 *(1-nRef_k4) - betaNref_k4 * nRef_k4);
        k4h = (alphaHref_k4 *(1-hRef_k4) - betaHref_k4 * hRef_k4);

        % Update values using Runge-Kutta 4th order method
        Vref(i+1) = Vref(i)*(1 + frozenNoiseRef(i)) + dt * (k1V/6 + k2V/3 + k3V/3 + k4V/6);
        mRef(i+1) = mRef(i)                         + dt * (k1m/6 + k2m/3 + k3m/3 + k4m/6);
        nRef(i+1) = nRef(i)                         + dt * (k1n/6 + k2n/3 + k3n/3 + k4n/6);
        hRef(i+1) = hRef(i)                         + dt * (k1h/6 + k2h/3 + k3h/3 + k4h/6);
    end

    % spike times
    thresh   = 0;    % volts
    distance = 50;   % samples
    [~,spikeTimesRef] = findpeaks(Vref,'MinPeakHeight',thresh,'MinPeakDistance',distance);
    
    %%
    
    if strcmp(inputType,'deltaFunc')
        waveform = zeros(72,1);
        waveform(30) = -3;
    elseif strcmp(inputType,'datasetEspike')
        waveform = [0
                    0.0091
                    0.0173
                    0.0252
                    0.0328
                    0.0397
                    0.0460
                    0.0515
                    0.0558
                    0.0590
                    0.0616
                    0.0637
                    0.0650
                    0.0651
                    0.0632
                    0.0593
                    0.0543
                    0.0491
                    0.0434
                    0.0358
                    0.0237
                    0.0077
                   -0.0064
                   -0.0237
                   -0.0492
                   -0.0658
                   -0.0803
                   -0.1124
                   -0.1355
                   -0.1612
                   -0.2359
                   -0.4428
                   -1.2874
                   -3.3540
                   -5.8737
                   -6.9707
                   -5.9490
                   -3.8268
                   -1.7193
                   -0.0983
                    0.8877
                    1.3374
                    1.4860
                    1.4537
                    1.2955
                    1.0996
                    0.9166
                    0.7494
                    0.6041
                    0.4871
                    0.3933
                    0.3157
                    0.2509
                    0.1963
                    0.1509
                    0.1143
                    0.0846
                    0.0599
                    0.0396
                    0.0225
                    0.0069
                   -0.0073
                   -0.0202
                   -0.0322
                   -0.0438
                   -0.0551
                   -0.0669
                   -0.0792
                   -0.0913
                   -0.1028
                   -0.1142
                   -0.1259];
    end
           
%     waveform = [0
%                -0.0038
%                -0.0142
%                -0.0073
%                 0.0040
%                 0.0074
%                 0.0186
%                 0.0312
%                 0.0277
%                 0.0149
%                 0.0031
%                -0.0021
%                 0.0052
%                 0.0108
%                -0.0014
%                -0.0156
%                -0.0182
%                -0.0209
%                -0.0310
%                -0.0389
%                -0.0372
%                -0.0380
%                -0.0583
%                -0.0843
%                -0.0862
%                -0.0663
%                -0.0518
%                -0.0533
%                -0.0613
%                -0.0721
%                -0.0867
%                -0.0910
%                -0.0802
%                -0.0918
%                -0.1471
%                -0.1833
%                -0.1206
%                 0.0227
%                 0.1153
%                 0.0134
%                -0.2700
%                -0.5157
%                -0.5345
%                -0.3842
%                -0.2337
%                -0.1413
%                -0.0866
%                -0.0733
%                -0.0979
%                -0.1169
%                -0.1056
%                -0.0899
%                -0.0892
%                -0.0796
%                -0.0561
%                -0.0567
%                -0.0783
%                -0.0801
%                -0.0666
%                -0.0611
%                -0.0554
%                -0.0499
%                -0.0559
%                -0.0558
%                -0.0397
%                -0.0280
%                -0.0238
%                -0.0140
%                 0.0035
%                 0.0277
%                 0.0485
%                 0.0456];

           %     waveform = [1; zeros(23,1)];
    waveform = [waveform; waveform(end)*((length(waveform):-1:1)'/length(waveform))];
    
    % 
    I_wave = ones(1,length(t))*0.0;
    
%     UL = length(t) - length(waveform) - 1;
%     LL = 30;
%     numElements = UL-LL+1;
%     r = randperm(numElements, nWaveTrials);
%     waveIdx = r + LL;
    
    waveIdx = spikeTimesRef;
    I_wave(waveIdx) = 1;

    I_wave = gain*conv(I_wave,waveform);
    I_wave = I_wave(1:length(t));
    
    I_current_tar = I_current_tar - I_wave;
    I_current_tar(1,50*(1/dt)) = 0; % reset first fifty ms to 0
    
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
    % Loop over the second set of neurons
    for i = 1:length(t)-1
        % K1
        [alphaM, betaM] = m_equations(V(i), Vrest);
        [alphaN, betaN] = n_equations(V(i), Vrest);
        [alphaH, betaH] = h_equations(V(i), Vrest);
        conductance_K = g_K*(n(i)^4);
        conductance_Na = g_Na*(m(i)^3)*h(i);
        I_Na(i) = conductance_Na*(V(i)-E_Na);
        I_K(i) = conductance_K*(V(i)-E_K);
        I_Leak(i) = g_Leak*(V(i)-E_Leak);
        Input = juncSwitch(i)*I_current_tar(i) - (I_Na(i) + I_K(i) + I_Leak(i));
        k1V = Input / C;
        k1m = (alphaM *(1-m(i)) - betaM * m(i));
        k1n = (alphaN *(1-n(i)) - betaN * n(i));
        k1h = (alphaH *(1-h(i)) - betaH * h(i));

        % K2
        Vk2 = V(i) + 0.5 * dt * k1V;
        mk2 = m(i) + 0.5 * dt * k1m;
        nk2 = n(i) + 0.5 * dt * k1n;
        hk2 = h(i) + 0.5 * dt * k1h;
        [alphaM_k2, betaM_k2] = m_equations(Vk2, Vrest);
        [alphaN_k2, betaN_k2] = n_equations(Vk2, Vrest);
        [alphaH_k2, betaH_k2] = h_equations(Vk2, Vrest);
        conductance_K_k2 = g_K*(nk2^4);
        conductance_Na_k2 = g_Na*(mk2^3)*hk2;
        I_Na_k2 = conductance_Na_k2*(Vk2-E_Na);
        I_K_k2 = conductance_K_k2*(Vk2-E_K);
        I_Leak_k2 = g_Leak*(Vk2-E_Leak);
        Input_k2 = I_current_tar(i) - (I_Na_k2 + I_K_k2 + I_Leak_k2);
        k2V = Input_k2 / C;
        k2m = (alphaM_k2 *(1-mk2) - betaM_k2 * mk2);
        k2n = (alphaN_k2 *(1-nk2) - betaN_k2 * nk2);
        k2h = (alphaH_k2 *(1-hk2) - betaH_k2 * hk2);

        % K3
        Vk3 = V(i) + 0.5 * dt * k2V;
        mk3 = m(i) + 0.5 * dt * k2m;
        nk3 = n(i) + 0.5 * dt * k2n;
        hk3 = h(i) + 0.5 * dt * k2h;
        [alphaM_k3, betaM_k3] = m_equations(Vk3, Vrest);
        [alphaN_k3, betaN_k3] = n_equations(Vk3, Vrest);
        [alphaH_k3, betaH_k3] = h_equations(Vk3, Vrest);
        conductance_K_k3 = g_K*(nk3^4);
        conductance_Na_k3 = g_Na*(mk3^3)*hk3;
        I_Na_k3 = conductance_Na_k3*(Vk3-E_Na);
        I_K_k3 = conductance_K_k3*(Vk3-E_K);
        I_Leak_k3 = g_Leak*(Vk3-E_Leak);
        Input_k3 = I_current_tar(i) - (I_Na_k3 + I_K_k3 + I_Leak_k3);
        k3V = Input_k3 / C;
        k3m = (alphaM_k3 *(1-mk3) - betaM_k3 * mk3);
        k3n = (alphaN_k3 *(1-nk3) - betaN_k3 * nk3);
        k3h = (alphaH_k3 *(1-hk3) - betaH_k3 * hk3);

        % K4
        Vk4 = V(i) + dt * k3V;
        mk4 = m(i) + dt * k3m;
        nk4 = n(i) + dt * k3n;
        hk4 = h(i) + dt * k3h;
        [alphaM_k4, betaM_k4] = m_equations(Vk4, Vrest);
        [alphaN_k4, betaN_k4] = n_equations(Vk4, Vrest);
        [alphaH_k4, betaH_k4] = h_equations(Vk4, Vrest);
        conductance_K_k4 = g_K*(nk4^4);
        conductance_Na_k4 = g_Na*(mk4^3)*hk4;
        I_Na_k4 = conductance_Na_k4*(Vk4-E_Na);
        I_K_k4 = conductance_K_k4*(Vk4-E_K);
        I_Leak_k4 = g_Leak*(Vk4-E_Leak);
        Input_k4 = I_current_tar(i) - (I_Na_k4 + I_K_k4 + I_Leak_k4);
        k4V = Input_k4 / C;
        k4m = (alphaM_k4 *(1-mk4) - betaM_k4 * mk4);
        k4n = (alphaN_k4 *(1-nk4) - betaN_k4 * nk4);
        k4h = (alphaH_k4 *(1-hk4) - betaH_k4 * hk4);

        % Update values using Runge-Kutta 4th order method
        V(i+1) = V(i)*(1 + frozenNoiseTar(i)) + dt * (k1V/6 + k2V/3 + k3V/3 + k4V/6);
        m(i+1) = m(i)                         + dt * (k1m/6 + k2m/3 + k3m/3 + k4m/6);
        n(i+1) = n(i)                         + dt * (k1n/6 + k2n/3 + k3n/3 + k4n/6);
        h(i+1) = h(i)                         + dt * (k1h/6 + k2h/3 + k3h/3 + k4h/6);
    end

    
    % spike times
    thresh   = 0;    % volts
    distance = 50;   % samples
    [~,spikeTimes] = findpeaks(V,'MinPeakHeight',thresh,'MinPeakDistance',distance);
%     
%     Vthresh = 15; 
%     Ithresh = 0;  
%     I = (-I_current_tar - (I_Na + I_K + I_Leak));
%     
%     % Find time points where voltage crosses the threshold during sodium upswing
%     above_threshold_indices = find((V(2:end) > Vthresh) & (I < Ithresh));
% 
%     % Identify sequential samples above threshold
%     diff_indices = diff(above_threshold_indices);
%     spikeTimes = [above_threshold_indices(1) above_threshold_indices(diff_indices > 5)];

    %%
    
    Duration = 0.010;
    
    stimTimes     = sort(waveIdx)'/fs;
    spikeTimesTar = spikeTimes'/fs;
    
    % bursting check
%     ISI    = diff(spikeTimesTar);
%     ISIbin = ISI < 0.008;
%     ISIbin = [false; ISIbin]; % correct for lost sample
% 
%     spikeTimesTarShortPreISI  = spikeTimesTar(ISIbin);
%     spikeTimesTarLongPreISI   = spikeTimesTar(~ISIbin);
    
    [ccg,tR] = CCG([stimTimes;spikeTimesTar],[ones(size(stimTimes));2*ones(size(spikeTimesTar))], ...
                    'binSize', 1/fs, 'duration', Duration, 'Fs', 1/fs, ...
                    'norm', 'counts');
    
     
    filterFlag     = false; 
    noDemeanFlag   = true;
    fpass          = 300;
    preLength      = 5*30;
    postLength     = 5*30;
    evokedSpikelet = waveformAvg(V',stimTimes(stimTimes < 60)*fs,preLength,postLength,fpass,fs,filterFlag,noDemeanFlag);            

%     figure
%     tiledlayout(2,1)

    if strcmp(inputType,'deltaFunc')
        nexttile(6)
    elseif strcmp(inputType,'datasetEspike')
        nexttile(5)
    end
    plot(tR*1e3,ccg(:,1,2)/length(spikeTimesRef),'k','LineWidth',2)
    xlim([0 2])
    xlabel('[ms]')
    ylabel('Prob.')
%     title('R v T')
    set(gca,'FontSize',6)
    box off
    
    if strcmp(inputType,'deltaFunc')
        nexttile(12)
    elseif strcmp(inputType,'datasetEspike')
        nexttile(11)
    end
    plot((-(preLength-1):postLength)*dt,evokedSpikelet*1000,'LineWidth',2) 
    xlim([0 2])
    xlabel('[ms]')
    ylabel('Vm [mV]')
% 	title('pulse')
    set(gca,'FontSize',6)
%     set(gca, 'YDir','reverse')
    box off
    
    if strcmp(inputType,'deltaFunc')
        nexttile(18)
    elseif strcmp(inputType,'datasetEspike')
        nexttile(17)
    end
%     yyaxis left
    plot((-(preLength-1):(postLength-1))*dt,33*diff(evokedSpikelet)*1000,'-.','LineWidth',2,'color',"#0072BD")
%     yyaxis right
%     plot((0:143)*dt*1000,waveform,'LineWidth',4) 
    xlim([0 2])
    xlabel('[ms]')
    ylabel('dVm/dt [mV/ms]')    
%     title('d(pulse)/dt')
    set(gca,'FontSize',6)
%     set(gca, 'YDir','reverse')
    box off
    
%     figure
%     tiledlayout(2,1)
%     nexttile; plot(tR*1000,ccg(:,1,2)); xlim([-1 5])
%     nexttile; plot((1:72*2)*dt,waveform); xlim([-1 5])
%     
%     figure
%     yyaxis left
%     plot(V(1:1e4))
%     hold on
%     plot(V(1:1e4))
%     hold off
% 
%     yyaxis right
%     plot(I_current_tar(1:1e4))
%     
% %     tiledlayout(2,2)            
% %     nexttile; plot(tR*1000,ccg(:,1,1))
% %     nexttile; plot(tR*1000,ccg(:,1,2))
% %     nexttile; plot(tR*1000,ccg(:,2,1))
% %     nexttile; plot(tR*1000,ccg(:,2,2))
%     
%     figure
%     yyaxis left
%     plot(V(1+1e6:1e4+1e6))
% 
%     yyaxis right
%     plot(I_current_tar(1+1e6:1e4+1e6))
    
    %%
    
%     tiledlayout(5,1)
%     
%     nexttile
%     plot(t,I_current_tar)
%     
%     nexttile
%     plot(t,V(:,1:end-1))
%     
%     binSize        = 1/fs;
%     trialLength    = 60; % ms
%     preTrialLength = 30; % ms
%     
%     spikeTimes     = spikeTimes*(1/30);
%     trialTimes     = [(waveIdx - preTrialLength)*(1/30);(waveIdx + trialLength)*(1/30)];
%     
%     color = 'k';
    
%     nexttile
%     [histCounts] = PSTH(spikeTimes',trialTimes',preTrialLength,trialLength,binSize,color);
%     plot((-29:1:60)/30,histCounts)
%     
%     nexttile
%     plot((0:23)/30,waveform)
%     xlim([-1 2])
%     
%     nexttile 
%     plot(t,I_Na+I_K+I_Leak+I_current)
    
    % removing the synchronous spikes                
%     [CCGtracked,spikeTimesInBin] = CCGtrackedSpikeTimes(stimTimes',spikeTimesTar',1/fs,Duration/2);
% 
%     refSpon = spikeTimesInBin{230}(:,1);
%     tarSpon = spikeTimesInBin{230}(:,2);
%     
%     refTrough = spikeTimesInBin{181}(:,1);
%     tarTrough = spikeTimesInBin{181}(:,2);
%     
%     refPeak = spikeTimesInBin{180}(:,1);
%     tarPeak = spikeTimesInBin{180}(:,2);
% 
%     filterFlag = false; 
%     fpass      = 300;
%     preLength  = 36;
%     postLength = 36;
%     
%     ACG      = CCGtrackedSpikeTimes(spikeTimesTar',spikeTimesTar',1/fs, Duration/2);
%     ACGsynch = CCGtrackedSpikeTimes(stimTimes',    tarPeak',      1/fs, Duration/2);
%     
%     ACG(601) = 0;
%     ACGsynch(601) = 0;
%     
%     plot(tR*1000,ACG)
%     hold on
%     plot(tR*1000,ACGsynch)
%     hold off
%     
% %     I = (I_current_tar - (I_Na + I_K + I_Leak))';
%     I = -I_Na';
%     I = (-(I_Na + I_K + I_Leak))';
%     
%     noDemeanFlag = true;
%     
%     [Ispon,  IsponAll]   = waveformAvg(I',tarSpon*fs,preLength,postLength,fpass,fs,filterFlag,noDemeanFlag);  
%     [Itrough,ItroughAll] = waveformAvg(I',tarTrough*fs,preLength,postLength,fpass,fs,filterFlag,noDemeanFlag);
%     [Ipeak,  IpeakAll]   = waveformAvg(I',tarPeak*fs,preLength,postLength,fpass,fs,filterFlag,noDemeanFlag);
%     
%     [Vspon,  VsponAll]   = waveformAvg(V',tarSpon*fs,preLength,postLength,fpass,fs,filterFlag,noDemeanFlag);  
%     [Vtrough,VtroughAll] = waveformAvg(V',tarTrough*fs,preLength,postLength,fpass,fs,filterFlag,noDemeanFlag);
%     [Vpeak,  VpeakAll]   = waveformAvg(V',tarPeak*fs,preLength,postLength,fpass,fs,filterFlag,noDemeanFlag);
%     
%     plot([Ispon Itrough Ipeak],'LineWidth',2)
%     legend('Ispon','Itrough','Ipeak')
%     
%     plot([Vspon Vtrough Vpeak],'LineWidth',2)
%     legend('Vspon','Vtrough','Vpeak')
%     
%     plot(Vspon,  Ispon,  'LineWidth',2)
%     hold on 
%     plot(Vtrough,Itrough,'LineWidth',2)
%     plot(Vpeak,  Ipeak,  'LineWidth',2)
%     hold off
%     legend('Ispon','Itrough','Ipeak')
    
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