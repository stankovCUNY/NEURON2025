function EIFephapse(timeMin,inputType)
    
    tic

%     close all
    
    time = timeMin*60;% in seconds
    
    noiseSc = 300; % normal noise scale

    Vrest = -65; % mV − change this to −65 if desired
    dt = 1/30; % ms
    fs  = (1/dt)*1000; % hz
    totalTime = time*1e3; % ms
    
    % Neuron parameters
    C = 281;                 % Membrane capacitance (pF)
    gL = 30;                 % Leak conductance (nS)
    EL = -60;              % Leak reversal potential (mV)
    DeltaT = 0.1;              % Slope factor (mV)
    VT = -50.4;              % Spike threshold (mV)
    V_reset = -65;         % Reset potential (mV)
    
    gain = 1500;

    % Initialization
    V = EL;                  % Membrane potential (mV)

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
    
    I_current_ref(50/dt:end) = 300; % Input of 3 microA/cm2 beginning at 50 ms and steady until end of time period.
    I_current_tar((2*60*1000)/dt:end) = 300; % Input of 3 microA/cm2 beginning at 50 ms and steady until end of time period.
    
    % add noise
    I_current_ref = I_current_ref + frozenNoiseRef;
    I_current_tar = I_current_tar + frozenNoiseTar;
    
    spikeTimesRef = [];
    for i = 1:length(t)

        dVdt = (-gL*(V - EL) + gL*DeltaT*exp((V - VT)/DeltaT) + I_current_ref(i))/C;
%         V = V*(1+frozenNoiseRef(i)) + dVdt * dt;
        V = V + dVdt*dt;

        if V >= VT
            spikeTimesRef = [spikeTimesRef i];
            V = V_reset;
        end

        V_trace(i) = V;
    end
    
    %%
    
    if strcmp(inputType,'deltaFunc')
        waveform = zeros(72,1);
        waveform(30) = -3;
    elseif strcmp(inputType,'datasetEspike')
        waveform = 0.1*[0
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
    
    spikeTimesTar = [];
    for i = 1:length(t)

        dVdt = (-gL*(V - EL) + gL*DeltaT*exp((V - VT)/DeltaT) + I_current_tar(i))/C;
%         V = V*(1+frozenNoiseTar(i)) + dVdt * dt;
        V = V + dVdt*dt;

        if V >= VT
            spikeTimesTar = [spikeTimesTar i];
            V = V_reset;
        end

        V_trace(i) = V;
    end
    
    %%
    
    Duration = 0.010;
    
    stimTimes     = sort(waveIdx)'/fs;
    spikeTimesTar = spikeTimesTar'/fs;
    
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
    evokedSpikelet = waveformAvg(V_trace',stimTimes(stimTimes < 60)*fs,preLength,postLength,fpass,fs,filterFlag,noDemeanFlag);            

%     figure
%     tiledlayout(2,1)

    if strcmp(inputType,'deltaFunc')
        nexttile(6)
    elseif strcmp(inputType,'datasetEspike')
        nexttile(5)
    end
    plot(tR*1e3,ccg(:,1,2)/length(spikeTimesRef),'k','LineWidth',1)
    xlim([0 2])
    xlabel('[ms]')
    ylabel('Prob.')
%     title('R v T')
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
% 	title('pulse')
    set(gca,'FontSize',5)
    set(gca,'FontName','Arial')
%     set(gca, 'YDir','reverse')
    box off
    
    if strcmp(inputType,'deltaFunc')
        nexttile(18)
    elseif strcmp(inputType,'datasetEspike')
        nexttile(17)
    end
%     yyaxis left
    plot((-(preLength-1):(postLength-1))*dt,33*diff(evokedSpikelet)*1000,'-.','LineWidth',1,'color',"#0072BD")
%     yyaxis right
%     plot((0:143)*dt*1000,waveform,'LineWidth',4) 
    xlim([0 2])
    xlabel('[ms]')
    ylabel('dVm/dt [mV/ms]')    
%     title('d(pulse)/dt')
    set(gca,'FontSize',5)
    set(gca,'FontName','Arial')
%     set(gca, 'YDir','reverse')
    box off
    
    %%
    
    figure(); histogram(V_trace(3e4*60*2:end))
    
    avgFRtar = length(spikeTimesTar)/(12*60)
     
%     binSize        = 1/fs;
%     lagLimit       = Duration/2;
%     
%     [CCGtracked,spikeTimesInBin] = CCGtrackedSpikeTimes(stimTimes,spikeTimesTar,binSize,lagLimit);
%     
%     Vm = [];
%     tRlist = [];
% 
%     for loopBin = 1:size(spikeTimesInBin,2)
% 
%        % convert from spike times to index
%        idx = round(spikeTimesInBin{1,loopBin}(:,1)*fs);
% 
%        Vm = [V(idx-0) Vm];
%        tRlist = [tR(loopBin)*ones(1,length(idx)) tRlist];
% 
%     end
% 
%     scatter(tRlist*1e3,Vm,'.')
%     ylim([-75 -50])
%     xlabel('t  [ms]')
%     ylabel('i Vm [mV]')

    
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
