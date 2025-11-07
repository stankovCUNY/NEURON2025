function HodgkinHuxleyEphapseCommonInput(timeMin,inputType)
    
    tic

%     close all
    
    time = timeMin*60;% in seconds
    
    noiseSc = 2.5; % normal noise scale

    Vrest = -65; % mV − change this to −65 if desired
    dt = 1/30; % ms
    fs  = (1/dt)*1000; % hz
    totalTime = time*1e3; % ms
    C = 1; % uF/cm^2
    
    % constants; values based on Table 1
%     E_Na = 115 + Vrest; % mV
%     E_K = -6 + Vrest; %mV
%     E_Leak = 10.6 + Vrest; % mV
%     
%     g_Na = 120; % mS/cm^2
%     g_K = 36; % mS/cm^2
%     g_Leak = 0.3; % mS/cm^2
    
    E_Na   = 55; % mV
    E_K    = -90; %mV
    E_Leak = -65; % mV
    
    g_Na = 35; % mS/cm^2
    g_K = 9; % mS/cm^2
    g_Leak = 0.1; % mS/cm^2

    gain = 6;
    
    % Vector of timesteps
    t = [0:dt:totalTime];
    
    % target ion channels switch
    ionSwitch = ones(1,length(t));
    ionSwitch(1:60*fs) = 0;
    
    % junction switch
    juncSwitch = ones(1,length(t));
    juncSwitch(60*fs+1:120*fs) = 0;
    
    frozenNoiseCI  = noiseSc*normrnd(0,1,1,length(t));
    frozenNoiseRef = noiseSc*normrnd(0,1,1,length(t));
    frozenNoiseTar = noiseSc*normrnd(0,1,1,length(t));
    
    %%
    
     % Current input −− change this to see how different inputs affect the neuron
    I_current_CI  = ones(1,length(t))*0.0; 
    I_current_ref = ones(1,length(t))*0.0;
    I_current_tar = ones(1,length(t))*0.0;
    
    I_current_CI(50/dt:end)           = 5;
    I_current_ref(50/dt:end)          = 5; % Input of 3 microA/cm2 beginning at 50 ms and steady until end of time period.
    I_current_tar((2*60*1000)/dt:end) = 5; % Input of 3 microA/cm2 beginning at 50 ms and steady until end of time period.
    
    % add noise
    I_current_CI  = I_current_CI  + frozenNoiseCI;
    I_current_ref = I_current_ref + frozenNoiseRef;
    I_current_tar = I_current_tar + frozenNoiseTar;
    
    % initializing values
    Vci(1)  = Vrest; % membrane potential is starting at its resting state
    Vref(1) = Vrest; 
    Vtar(1) = Vrest;
    s(1)    = 0; % initial synaptic gain
    
    % separate functions to get the alpha and beta values
    [alphaMci, betaMci]   = m_equations(Vci(1), Vci);
    [alphaNci, betaNci]   = n_equations(Vci(1), Vci);
    [alphaHci, betaHci]   = h_equations(Vci(1), Vci);
    
    [alphaMref, betaMref] = m_equations(Vref(1), Vrest);
    [alphaNref, betaNref] = n_equations(Vref(1), Vrest);
    [alphaHref, betaHref] = h_equations(Vref(1), Vrest);
    
    [alphaMtar, betaMtar] = m_equations(Vtar(1), Vrest);
    [alphaNtar, betaNtar] = n_equations(Vtar(1), Vrest);
    [alphaHtar, betaHtar] = h_equations(Vtar(1), Vrest);
    
    % initializing gating variables to the asymptotic values when membrane potential
    % is set to the membrane resting value based on equation 13
    mCI(1)  = (alphaMci / (alphaMci + betaMci));
    nCI(1)  = (alphaNci / (alphaNci + betaNci));
    hCI(1)  = (alphaHci / (alphaHci + betaHci));
    
    mRef(1) = (alphaMref / (alphaMref + betaMref));
    nRef(1) = (alphaNref / (alphaNref + betaNref));
    hRef(1) = (alphaHref / (alphaHref + betaHref));
    
    mTar(1) = (alphaMtar / (alphaMtar + betaMtar));
    nTar(1) = (alphaNtar / (alphaNtar + betaNtar));
    hTar(1) = (alphaHtar / (alphaHtar + betaHtar));
    
    % solve for common input spike times
    for i = 1:length(t)
        
        % calculate new alpha and beta based on last known membrane potenatial
        [alphaMci, betaMci]   = m_equations(Vci(i), Vrest);
        [alphaNci, betaNci]   = n_equations(Vci(i), Vrest);
        [alphaHci, betaHci]   = h_equations(Vci(i), Vrest);
        
        % conductance variables − computed separately to show how this
        % changes with membrane potential in one of the graphs
        conductance_K_ci(i)   = g_K*(nCI(i)^4);
        conductance_Na_ci(i)  = g_Na*(mCI(i)^3)*hCI(i);
        
        % retrieving ionic currents
        I_Na_ci(i)   = conductance_Na_ci(i)*(Vci(i)-E_Na);
        I_K_ci(i)    = conductance_K_ci(i)*(Vci(i)-E_K);
        I_Leak_ci(i) = g_Leak*(Vci(i)-E_Leak);
        
        % Calculating the input
%         Input_ci  = - (I_Na_ci(i) + I_K_ci(i) + I_Leak_ci(i));
        Input_ci  = I_current_CI(i) - (I_Na_ci(i) + I_K_ci(i) + I_Leak_ci(i));
        
        % Calculating the new membrane potential
%         Vci(i+1)  = Vci(i)*(1 + frozenNoiseCI(i))  + Input_ci*  dt*(1/C); 
        Vci(i+1)  = Vci(i)  + Input_ci*dt*(1/C); 

        % getting new values for the gating variables
        mCI(i+1)  = mCI(i) + (alphaMci *(1-mCI(i)) - betaMci * mCI(i))*dt;
        nCI(i+1)  = nCI(i) + (alphaNci *(1-nCI(i)) - betaNci * nCI(i))*dt;
        hCI(i+1)  = hCI(i) + (alphaHci *(1-hCI(i)) - betaHci * hCI(i))*dt;
        
    end
    
    % spike times
    thresh   = 0;    % volts
    distance = 10;   % samples
    [~,spikeTimesCI] = findpeaks(Vci,'MinPeakHeight',thresh,'MinPeakDistance',distance);
    
    % ramdom stims
%     desiredStimFreq = 500;
%     numStims = desiredStimFreq*((timeMin*60));
%     spikeTimesCI = sort(randperm(length(t),numStims));

    %%
    
%     waveform = zeros(72,1); 
%     waveform(30:32) = -3;
%     waveform(32:52) = 0.3;
    
%     waveform = [zeros(36,1); -0.05*exp(-(0:35)*0.1)'];
    
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
    
    waveIdx = spikeTimesCI;
    I_wave(waveIdx) = 1;

    I_wave = gain*conv(I_wave,waveform);
    I_wave = I_wave(1:length(t));
    
    % Current input −− change this to see how different inputs affect the neuron
%     I_current_ref = ones(1,length(t))*0.0;
%     I_current_tar = ones(1,length(t))*0.0;
    
    % Current input with sine wave
%     Fc = 8;
%     I_current_ref = 7.5*sin(2*pi*Fc*t/1000); % t is in ms so divide by 1000
%     I_current_tar = 7.5*sin(2*pi*Fc*t/1000); % t is in ms so divide by 1000
    
%     I_current_ref(50/dt:end) = 0; % Input of 1 microA/cm2 beginning at 50 ms and steady until end of time period.
%     I_current_tar(30/dt:end) = 0; % Input of 1 microA/cm2 beginning at 30 ms and steady until end of time period.
    
    I_current_ref = I_current_ref - I_wave;
%     I_current_ref(1,50*(1/fs)) = 0; % reset first fifty ms to 0
    
    I_current_tar = I_current_tar - I_wave;
%     I_current_tar(1,50*(1/fs)) = 0; % reset first fifty ms to 0
    
    % repeat for time determined in totalTime , by each dt
    for i = 1:length(t)
        
        % calculate new alpha and beta based on last known membrane potenatial
        [alphaMref, betaMref] = m_equations(Vref(i), Vrest);
        [alphaNref, betaNref] = n_equations(Vref(i), Vrest);
        [alphaHref, betaHref] = h_equations(Vref(i), Vrest);
    
        [alphaMtar, betaMtar] = m_equations(Vtar(i), Vrest);
        [alphaNtar, betaNtar] = n_equations(Vtar(i), Vrest);
        [alphaHtar, betaHtar] = h_equations(Vtar(i), Vrest);
        
        % conductance variables − computed separately to show how this
        % changes with membrane potential in one of the graphs
        conductance_K_ref(i)  = g_K*(nRef(i)^4);
        conductance_Na_ref(i) = g_Na*(mRef(i)^3)*hRef(i);
        
        conductance_K_tar(i)  = g_K*(nTar(i)^4);
        conductance_Na_tar(i) = g_Na*(mTar(i)^3)*hTar(i);
        
        % retrieving ionic currents
        I_Na_ref(i)   = conductance_Na_ref(i)*(Vref(i)-E_Na);
        I_K_ref(i)    = conductance_K_ref(i)*(Vref(i)-E_K);
        I_Leak_ref(i) = g_Leak*(Vref(i)-E_Leak);
        
        I_Na_tar(i)   = conductance_Na_tar(i)*(Vtar(i)-E_Na);
        I_K_tar(i)    = conductance_K_tar(i)*(Vtar(i)-E_K);
        I_Leak_tar(i) = g_Leak*(Vtar(i)-E_Leak);
        
        % Calculating the input
        Input_ref = juncSwitch(i)*(I_current_ref(i)) - (ionSwitch(i)*(I_Na_ref(i) + I_K_ref(i)) + I_Leak_ref(i));
        Input_tar = juncSwitch(i)*(I_current_tar(i)) - (ionSwitch(i)*(I_Na_tar(i) + I_K_tar(i)) + I_Leak_tar(i));
        
        % Calculating the new membrane potential
%         Vref(i+1) = Vref(i)*(1 + frozenNoiseRef(i)) + Input_ref* dt*(1/C);
%         Vtar(i+1) = Vtar(i)*(1 + frozenNoiseTar(i)) + Input_tar* dt*(1/C);
        Vref(i+1) = Vref(i) + Input_ref*dt*(1/C);
        Vtar(i+1) = Vtar(i) + Input_tar*dt*(1/C);
        
        % getting new values for the gating variables
        mRef(i+1) = mRef(i) + (alphaMref *(1-mRef(i)) - betaMref * mRef(i))*dt;
        nRef(i+1) = nRef(i) + (alphaNref *(1-nRef(i)) - betaNref * nRef(i))*dt;
        hRef(i+1) = hRef(i) + (alphaHref *(1-hRef(i)) - betaHref * hRef(i))*dt;
        
        mTar(i+1) = mTar(i) + (alphaMtar *(1-mTar(i)) - betaMtar * mTar(i))*dt;
        nTar(i+1) = nTar(i) + (alphaNtar *(1-nTar(i)) - betaNtar * nTar(i))*dt;
        hTar(i+1) = hTar(i) + (alphaHtar *(1-hTar(i)) - betaHtar * hTar(i))*dt;
        
    end
    
    % spike times
    thresh   = 0;    % volts
    distance = 50;   % samples
    [~,spikeTimesRef] = findpeaks(Vref,'MinPeakHeight',thresh,'MinPeakDistance',distance);
    [~,spikeTimesTar] = findpeaks(Vtar,'MinPeakHeight',thresh,'MinPeakDistance',distance);
    
%     Vthresh = 15; 
%     Ithresh = 0;  
    
%     %% spike times ref
%     I = juncSwitch.*(I_current_ref) - (ionSwitch.*(I_Na_ref + I_K_ref) + I_Leak_ref);
%     V = Vref;
%     
%     % Find time points where voltage crosses the threshold during sodium upswing
%     above_threshold_indices = find((V(2:end) > Vthresh) & (I > Ithresh));
% 
%     % Identify sequential samples above threshold
%     diff_indices = diff(above_threshold_indices);
%     spikeTimesRef = [above_threshold_indices(1) above_threshold_indices(diff_indices > 5)];
%     
%     %% spike times tar
%     I = juncSwitch.*(I_current_tar) - (ionSwitch.*(I_Na_tar + I_K_tar) + I_Leak_tar);
%     V = Vtar;
%     
%     % Find time points where voltage crosses the threshold during sodium upswing
%     above_threshold_indices = find((V(2:end) > Vthresh) & (I > Ithresh));
% 
%     % Identify sequential samples above threshold
%     diff_indices = diff(above_threshold_indices);
%     spikeTimesTar = [above_threshold_indices(1) above_threshold_indices(diff_indices > 5)];
    
    %%
    
    Duration = 0.100;
    
    spikeTimesCI  = spikeTimesCI'/fs;
    spikeTimesRef = spikeTimesRef'/fs;
    spikeTimesTar = spikeTimesTar'/fs;
    
    [ccgRT,tR] = CCG([spikeTimesRef;spikeTimesTar],[ones(size(spikeTimesRef));2*ones(size(spikeTimesTar))], ...
                    'binSize', 1/fs, 'duration', Duration, 'Fs', 1/fs, ...
                    'norm', 'counts');
                
%     figure
%     tiledlayout(2,2)            
%     nexttile; plot(tR*1000,ccgRT(:,1,1))
%     title("R v R")
%     nexttile; plot(tR*1000,ccgRT(:,1,2))
%     title("R v T")
%     nexttile; plot(tR*1000,ccgRT(:,2,1))
%     title("T v R")
%     nexttile; plot(tR*1000,ccgRT(:,2,2))
%     title("T v T")
    
    [ccgCiR,tR] = CCG([spikeTimesCI;spikeTimesRef],[ones(size(spikeTimesCI));2*ones(size(spikeTimesRef))], ...
                    'binSize', 1/fs, 'duration', Duration, 'Fs', 1/fs, ...
                    'norm', 'counts');
                
%     figure
%     tiledlayout(2,2)            
%     nexttile; plot(tR*1000,ccgCiR(:,1,1))
%     title("CI v CI")
%     nexttile; plot(tR*1000,ccgCiR(:,1,2))
%     title("CI v R")
%     nexttile; plot(tR*1000,ccgCiR(:,2,1))
%     title("R v CI")
%     nexttile; plot(tR*1000,ccgCiR(:,2,2))
%     title("R v R")
    
    [ccgCiT,tR] = CCG([spikeTimesCI;spikeTimesTar],[ones(size(spikeTimesCI));2*ones(size(spikeTimesTar))], ...
                    'binSize', 1/fs, 'duration', Duration, 'Fs', 1/fs, ...
                    'norm', 'counts');
                
%     figure
%     tiledlayout(2,2)            
%     nexttile; plot(tR*1000,ccgCiT(:,1,1))
%     title('CI v CI')
%     nexttile; plot(tR*1000,ccgCiT(:,1,2))
%     title("CI v T")
%     nexttile; plot(tR*1000,ccgCiT(:,2,1))
%     title("T v CI")
%     nexttile; plot(tR*1000,ccgCiT(:,2,2))
%     title("T v T")
%     
%     figure
%     yyaxis left
%     plot(Vref(1:1e4))
%     hold on
%     plot(Vtar(1:1e4))
%     hold off
% 
%     yyaxis right
%     plot(I_CI_ref(1:1e4))
%     hold on
%     plot(I_CI_tar(1:1e4))
%     hold off
%     
%     title("left ref/tar Vm, right ref/tar EPSP")
    
%     tiledlayout(1,2,'Padding','none','TileSpacing','compact')
    
    nexttile 
    plot(tR*1000,ccgCiT(:,1,2)/length(spikeTimesCI),'LineWidth',1)
    hold on
    plot(tR*1000,ccgCiR(:,1,2)/length(spikeTimesCI),'LineWidth',1)
    hold off
    xlim([-10 10])
    xlabel('[ms]')
    ylabel('Spike Probability')
%     title('CI v T and CI v R')
    set(gca,'FontSize',5)
    set(gca,'FontName','Arial')
    box off
    
    nexttile
    plot(tR*1000,ccgRT(:,1,2)/length(spikeTimesRef),'k','LineWidth',1)
    xlim([-10 10])
    xlabel('[ms]')
    ylabel('Spike Probability')
%     title('R v T')
    set(gca,'FontSize',5)
    set(gca,'FontName','Arial')
    box off
    
    %% inset - location of subpart on figure
%     xstart = 0.40;
%     xend   = 0.475;
%     ystart = 0.12;
%     yend   = 0.15;
%     axes('position',[xstart ystart xend-xstart yend-ystart ])
%     
%     plot(tR*1000,ccgCiT(:,1,2),'LineWidth',2)
%     hold on
%     plot(tR*1000,ccgCiR(:,1,2),'LineWidth',2)
%     hold off
%     xlim([0.5 1.5])
%     xlabel('[ms]')
%     ylabel('Counts')
% %     title('R v T')
%     set(gca,'FontSize',6)
%     box off
    
    toc
    
    %% raster analysis
    
%     spikeTimesCI  = spikeTimesCI( spikeTimesCI > 120);
%     spikeTimesRef = spikeTimesRef(spikeTimesRef > 120);
%     spikeTimesTar = spikeTimesTar(spikeTimesTar > 120);
%         
%     raster(spikeTimesCI,spikeTimesRef,spikeTimesTar,20)
    %% synch spike analysis
    
%     binSize  = 1/fs;
%     lagLimit = Duration/2;
%     
%     [~,spikeTimesInBin] = CCGtrackedSpikeTimes(spikeTimesCI',spikeTimesRef',binSize,lagLimit);
% 
% %     binOfInterest1 = 1536;
% %     binOfInterest2 = 1568;
%     binOfInterest1 = 1531;
%     binOfInterest2 = 1529;
%     
%     spikeTimesCIsync1  = spikeTimesInBin{binOfInterest1}(:,1);
%     spikeTimesRefsync1 = spikeTimesInBin{binOfInterest1}(:,2);
%     
%     spikeTimesCIsync2  = spikeTimesInBin{binOfInterest2}(:,1);
%     spikeTimesRefsync2 = spikeTimesInBin{binOfInterest2}(:,2);
%     
% %     stimPeak1idx = round(spikeTimesCIsync1*fs + 36);
% %     stimPeak2idx = round(spikeTimesCIsync2*fs + 36);
% 
%     stimPeak1idx = round(spikeTimesRefsync1*fs + 36);
%     stimPeak2idx = round(spikeTimesRefsync2*fs + 36);
%     
% %     VrefPeak1               = Vref(stimPeak1idx);
% %     conductance_Na_refPeak1 = conductance_Na_ref(stimPeak1idx);
% %     conductance_K_refPeak1  = conductance_K_ref(stimPeak1idx);
% %     
% %     VrefPeak2               = Vref(stimPeak2idx);
% %     conductance_Na_refPeak2 = conductance_Na_ref(stimPeak2idx);
% %     conductance_K_refPeak2  = conductance_K_ref(stimPeak2idx);
%     
% %     tiledlayout(3,1)
% %     
% %     nexttile
% %     histogram(VrefPeak1)
% %     hold on
% %     histogram(VrefPeak2)
% %     hold off
% %     
% %     nexttile
% %     histogram(conductance_Na_refPeak1)
% %     hold on
% %     histogram(conductance_Na_refPeak2)
% %     hold off
% %     
% %     nexttile
% %     histogram(conductance_K_refPeak1)
% %     hold on
% %     histogram(conductance_K_refPeak2)
% %     hold off
%     
%     filterFlag     = false; 
%     noDemeanFlag   = true;
%     fpass          = 300;
%     preLength      = 5*30;
%     postLength     = 5*30;
%     
%     [evokedSpikeletMean1,evokedSpikelet1] = waveformAvg(Vref',stimPeak1idx,preLength,postLength,fpass,fs,filterFlag,noDemeanFlag);            
%     [evokedSpikeletMean2,evokedSpikelet2] = waveformAvg(Vref',stimPeak2idx,preLength,postLength,fpass,fs,filterFlag,noDemeanFlag);  
%     
%     figure()
%     plot(evokedSpikelet1*1000)
%     figure()
%     plot(evokedSpikelet2*1000)
%     
%     figure()
%     plot(evokedSpikeletMean1*1000)
%     hold on 
%     plot(evokedSpikeletMean2*1000)
%     hold off
%     
% %     % remove sync spikes for waveforms
% %     spikeTimesCInosync  = setdiff(spikeTimesCI, spikeTimesCIsync);
% %     spikeTimesRefnosync = setdiff(spikeTimesRef,spikeTimesRefsync);
% %     
% %     [ccgCiRsynch,tR] = CCG([spikeTimesCI;spikeTimesRefsync],[ones(size(spikeTimesCI));2*ones(size(spikeTimesRefsync))], ...
% %                        'binSize', 1/fs, 'duration', Duration, 'Fs', 1/fs, ...
% %                        'norm', 'counts');
% %     
% %     figure()
% %     plot(tR*1000,ccgCiRsynch(:,1,2),'LineWidth',2)
                   
end

%% original

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

function raster(ref,tarOn,tarOff,lgmx)
    
    % convert to ms
    ref    = ref*1000;
    tarOn  = tarOn*1000;
    tarOff = tarOff*1000;
    
    xR = []; yR = []; xTon = []; yTon = [];  xToff = []; yToff = [];  
    for k = 1:numel(ref)
        rtem    = ref(ref>=ref(k)-lgmx & ref<=ref(k)+lgmx) - ref(k);
        ttemOn  = tarOn(tarOn >=ref(k)-lgmx & tarOn <=ref(k)+lgmx) - ref(k);
        ttemOff = tarOff(tarOff>=ref(k)-lgmx & tarOff<=ref(k)+lgmx) - ref(k);
        xR      = [xR;rtem];
        yR      = [yR;k.*ones(size(rtem))];
        xTon    = [xTon;ttemOn];
        yTon    = [yTon;k.*ones(size(ttemOn))];
        xToff   = [xToff;ttemOff];
        yToff   = [yToff;k.*ones(size(ttemOff))];
    end
    
    scatter(xToff,yToff,1,'filled','markerfacecolor','k'); hold on;
    hold on
    scatter(xTon,yTon,1,'filled','markerfacecolor','r'); hold on;
    hold off
    xlabel('Lag times [sec]')
    ylabel('Trial')
    xlim([-lgmx,lgmx])
    title('RASTER: Target spikes aligned with reference spikes')

end