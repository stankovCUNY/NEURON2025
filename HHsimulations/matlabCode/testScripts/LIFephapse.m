function LIFephapse(timeMin,inputType)
    
    tic
    
    time = timeMin*60;% in seconds
    
    noiseSc = 0.02; % normal noise scale

    Vrest = -65; % mV − change this to −65 if desired
    dt = 1/30; % ms
    fs  = (1/dt)*1000; % hz
    totalTime = time*1e3; % ms
    C = 1; % uF/cm^2
    
    % constants; values based on Table 1
    E_Leak = 10.6 + Vrest; % mV
    
    g_Leak = 0.3; % mS/cm^2
    
    gain = 6;
    
    threshold_value = -50;
    reset_value = -65; 
    
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
    
    spikeTimesRef = [];
    
    % Current input −− change this to see how different inputs affect the neuron
    I_current_ref = ones(1,length(t))*0.0;
    I_current_tar = ones(1,length(t))*0.0;
    
    I_current_ref(50/dt:end) = 0; % Input of 3 microA/cm2 beginning at 50 ms and steady until end of time period.
    I_current_tar(30/dt:end) = 0; % Input of 3 microA/cm2 beginning at 50 ms and steady until end of time period.
    
    % initializing values
    Vref(1)  = Vrest; % membrane potential is starting at its resting state
    
    % solve for common input spike times
    for i = 1:length(t)
        
        % Calculating the input
        Input_ref  = I_current_ref(i) - g_Leak*(Vref(i)-E_Leak);
        
        % Calculating the new membrane potential
        Vref(i+1)  = Vref(i)*(1 + frozenNoiseRef(i))  + Input_ref*  dt*(1/C); 
        
        % Check for spiking condition (threshold)
        if Vref(i+1) >= threshold_value
            % Spike occurred, reset membrane potential
            Vref(i+1) = reset_value;

            % Store spike time
            spikeTimesRef = [spikeTimesRef, t(i)];
        end
        
    end
    
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
           
    waveform = [waveform; waveform(end)*((length(waveform):-1:1)'/length(waveform))];
    
    I_wave = ones(1,length(t))*0.0;
    
    waveIdx = spikeTimesRef;
    I_wave(waveIdx) = 1;

    I_wave = gain*conv(I_wave,waveform);
    I_wave = I_wave(1:length(t));
    
    I_current_tar = I_current_tar - I_wave;
    I_current_tar(1,50*(1/dt)) = 0; % reset first fifty ms to 0
    
    % initializing values
    V(1) = Vrest; % membrane potential is starting at its resting state
    
    % repeat for time determined in totalTime , by each dt
    for i = 1:length(t)
        
        % Calculating the input
        Input = juncSwitch(i)*I_current_tar(i) - g_Leak*(V(i)-E_Leak);
        
        % Calculating the new membrane potential
        V(i+1) = V(i)*(1 + frozenNoiseTar(i)) + Input*dt*(1/C);
        
    end
    
    % spike times
    thresh   = 0;    % volts
    distance = 50;   % samples
    [~,spikeTimes] = findpeaks(V,'MinPeakHeight',thresh,'MinPeakDistance',distance);

    %%
    
    Duration = 0.010;
    
    stimTimes     = sort(waveIdx)'/fs;
    spikeTimesTar = spikeTimes'/fs;
        
    [ccg,tR] = CCG([stimTimes;spikeTimesTar],[ones(size(stimTimes));2*ones(size(spikeTimesTar))], ...
                    'binSize', 1/fs, 'duration', Duration, 'Fs', 1/fs, ...
                    'norm', 'counts');
    
     
    filterFlag     = false; 
    noDemeanFlag   = true;
    fpass          = 300;
    preLength      = 5*30;
    postLength     = 5*30;
    evokedSpikelet = waveformAvg(V',stimTimes(stimTimes < 60)*fs,preLength,postLength,fpass,fs,filterFlag,noDemeanFlag);            

    if strcmp(inputType,'deltaFunc')
        nexttile(6)
    elseif strcmp(inputType,'datasetEspike')
        nexttile(5)
    end
    plot(tR*1e3,ccg(:,1,2)/length(spikeTimesRef),'k','LineWidth',1)
    xlim([0 2])
    xlabel('[ms]')
    ylabel('Prob.')
    set(gca,'FontSize',5)
    set(gca,'FontName','Arial')
    box off
    
    if strcmp(inputType,'deltaFunc')
        nexttile(12)
    elseif strcmp(inputType,'datasetEspike')
        nexttile(11)
    end
    plot((-(preLength-1):postLength)*dt,evokedSpikelet*1000,'LineWidth',1) 
    xlim([0 2])
    xlabel('[ms]')
    ylabel('Vm [mV]')
    set(gca,'FontSize',5)
    set(gca,'FontName','Arial')
    box off
    
    if strcmp(inputType,'deltaFunc')
        nexttile(18)
    elseif strcmp(inputType,'datasetEspike')
        nexttile(17)
    end
    plot((-(preLength-1):(postLength-1))*dt,33*diff(evokedSpikelet)*1000,'-.','LineWidth',1,'color',"#0072BD")
    xlim([0 2])
    xlabel('[ms]')
    ylabel('dVm/dt [mV/ms]')    
    set(gca,'FontSize',5)
    set(gca,'FontName','Arial')
    box off
    
end
