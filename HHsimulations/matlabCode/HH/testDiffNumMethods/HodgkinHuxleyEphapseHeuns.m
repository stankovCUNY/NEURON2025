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
    
    % Heun's method integration loop for the reference neuron
    for i = 1:length(t)-1
        % Predictor step
        [alphaMref, betaMref] = m_equations(Vref(i), Vrest);
        [alphaNref, betaNref] = n_equations(Vref(i), Vrest);
        [alphaHref, betaHref] = h_equations(Vref(i), Vrest);
        conductance_K_ref = g_K*(nRef(i)^4);
        conductance_Na_ref = g_Na*(mRef(i)^3)*hRef(i);
        I_Na_ref(i) = conductance_Na_ref*(Vref(i)-E_Na);
        I_K_ref(i) = conductance_K_ref*(Vref(i)-E_K);
        I_Leak_ref(i) = g_Leak*(Vref(i)-E_Leak);
        Input_ref = I_current_ref(i) - (I_Na_ref(i) + I_K_ref(i) + I_Leak_ref(i));
        k1V = Input_ref / C + frozenNoiseRef(i);
        k1m = (alphaMref *(1-mRef(i)) - betaMref * mRef(i));
        k1n = (alphaNref *(1-nRef(i)) - betaNref * nRef(i));
        k1h = (alphaHref *(1-hRef(i)) - betaHref * hRef(i));

        Vref_predictor = Vref(i) + dt * k1V;
        mRef_predictor = mRef(i) + dt * k1m;
        nRef_predictor = nRef(i) + dt * k1n;
        hRef_predictor = hRef(i) + dt * k1h;

        % Corrector step
        [alphaMref_pred, betaMref_pred] = m_equations(Vref_predictor, Vrest);
        [alphaNref_pred, betaNref_pred] = n_equations(Vref_predictor, Vrest);
        [alphaHref_pred, betaHref_pred] = h_equations(Vref_predictor, Vrest);
        conductance_K_ref_pred = g_K*(nRef_predictor^4);
        conductance_Na_ref_pred = g_Na*(mRef_predictor^3)*hRef_predictor;
        I_Na_ref_pred = conductance_Na_ref_pred*(Vref_predictor-E_Na);
        I_K_ref_pred = conductance_K_ref_pred*(Vref_predictor-E_K);
        I_Leak_ref_pred = g_Leak*(Vref_predictor-E_Leak);
        Input_ref_pred = I_current_ref(i) - (I_Na_ref_pred + I_K_ref_pred + I_Leak_ref_pred);
        k2V = Input_ref_pred / C + frozenNoiseRef(i+1); % Use the predicted value for noise
        k2m = (alphaMref_pred *(1-mRef_predictor) - betaMref_pred * mRef_predictor);
        k2n = (alphaNref_pred *(1-nRef_predictor) - betaNref_pred * nRef_predictor);
        k2h = (alphaHref_pred *(1-hRef_predictor) - betaHref_pred * hRef_predictor);

        % Update values using Heun's method
        Vref(i+1) = Vref(i)*(1 + frozenNoiseRef(i)) + 0.5 * dt * (k1V + k2V);
        mRef(i+1) = mRef(i)                         + 0.5 * dt * (k1m + k2m);
        nRef(i+1) = nRef(i)                         + 0.5 * dt * (k1n + k2n);
        hRef(i+1) = hRef(i)                         + 0.5 * dt * (k1h + k2h);
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
    Vtar(1) = Vrest; % membrane potential is starting at its resting state
    
    % separate functions to get the alpha and beta values
    [alphaMtar, betaMtar] = m_equations(Vtar(1), Vrest);
    [alphaNtar, betaNtar] = n_equations(Vtar(1), Vrest);
    [alphaHtar, betaHtar] = h_equations(Vtar(1), Vrest);
    
    % initializing gating variables to the asymptotic values when membrane potential
    % is set to the membrane resting value based on equation 13
    mTar(1) = (alphaMtar / (alphaMtar + betaMtar));
    nTar(1) = (alphaNtar / (alphaNtar + betaNtar));
    hTar(1) = (alphaHtar / (alphaHtar + betaHtar));
    
    % repeat for time determined in totalTime , by each dt
    % Loop over the second set of neurons
    % Heun's method integration loop for the target neuron
    for i = 1:length(t)-1
        % Predictor step
        [alphaMtar, betaMtar] = m_equations(Vtar(i), Vrest);
        [alphaNtar, betaNtar] = n_equations(Vtar(i), Vrest);
        [alphaHtar, betaHtar] = h_equations(Vtar(i), Vrest);
        conductance_K_tar = g_K*(nTar(i)^4);
        conductance_Na_tar = g_Na*(mTar(i)^3)*hTar(i);
        I_Na_tar(i) = conductance_Na_tar*(Vtar(i)-E_Na);
        I_K_tar(i) = conductance_K_tar*(Vtar(i)-E_K);
        I_Leak_tar(i) = g_Leak*(Vtar(i)-E_Leak);
        Input_tar = I_current_tar(i) - (I_Na_tar(i) + I_K_tar(i) + I_Leak_tar(i));
        k1V = Input_tar / C + frozenNoiseTar(i);
        k1m = (alphaMtar *(1-mTar(i)) - betaMtar * mTar(i));
        k1n = (alphaNtar *(1-nTar(i)) - betaNtar * nTar(i));
        k1h = (alphaHtar *(1-hTar(i)) - betaHtar * hTar(i));

        Vtar_predictor = Vtar(i) + dt * k1V;
        mTar_predictor = mTar(i) + dt * k1m;
        nTar_predictor = nTar(i) + dt * k1n;
        hTar_predictor = hTar(i) + dt * k1h;

        % Corrector step
        [alphaMtar_pred, betaMtar_pred] = m_equations(Vtar_predictor, Vrest);
        [alphaNtar_pred, betaNtar_pred] = n_equations(Vtar_predictor, Vrest);
        [alphaHtar_pred, betaHtar_pred] = h_equations(Vtar_predictor, Vrest);
        conductance_K_tar_pred = g_K*(nTar_predictor^4);
        conductance_Na_tar_pred = g_Na*(mTar_predictor^3)*hTar_predictor;
        I_Na_tar_pred = conductance_Na_tar_pred*(Vtar_predictor-E_Na);
        I_K_tar_pred = conductance_K_tar_pred*(Vtar_predictor-E_K);
        I_Leak_tar_pred = g_Leak*(Vtar_predictor-E_Leak);
        Input_tar_pred = I_current_tar(i) - (I_Na_tar_pred + I_K_tar_pred + I_Leak_tar_pred);
        k2V = Input_tar_pred / C + frozenNoiseTar(i+1); % Use the predicted value for noise
        k2m = (alphaMtar_pred *(1-mTar_predictor) - betaMtar_pred * mTar_predictor);
        k2n = (alphaNtar_pred *(1-nTar_predictor) - betaNtar_pred * nTar_predictor);
        k2h = (alphaHtar_pred *(1-hTar_predictor) - betaHtar_pred * hTar_predictor);

        % Update values using Heun's method
        Vtar(i+1) = Vtar(i)*(1 + frozenNoiseTar(i)) + 0.5 * dt * (k1V + k2V);
        mTar(i+1) = mTar(i)                         + 0.5 * dt * (k1m + k2m);
        nTar(i+1) = nTar(i)                         + 0.5 * dt * (k1n + k2n);
        hTar(i+1) = hTar(i)                         + 0.5 * dt * (k1h + k2h);
    end

    % spike times
    thresh   = 0;    % volts
    distance = 50;   % samples
    [~,spikeTimes] = findpeaks(Vtar,'MinPeakHeight',thresh,'MinPeakDistance',distance);
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
    evokedSpikelet = waveformAvg(Vtar',stimTimes(stimTimes < 60)*fs,preLength,postLength,fpass,fs,filterFlag,noDemeanFlag);            

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