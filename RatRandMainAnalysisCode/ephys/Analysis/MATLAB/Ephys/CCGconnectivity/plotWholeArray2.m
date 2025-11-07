load('/home/nasko/CUNY_Work_Han_Kamran_Nat/data/wake_new/unitsGlobalSnapshots.mat')

channelLayout = [1 8 9  16 17 24 25 32 33 40 41 48 49 56 57 64 ...
                 2 7 10 15 18 23 26 31 34 39 42 47 50 55 58 63 ...
                 3 6 11 14 19 22 27 30 35 38 43 46 51 54 59 62 ...
                 4 5 12 13 20 21 28 29 36 37 44 45 52 53 60 62];

channelShankLayout = [  ones(8,1); 2*ones(8,1); 3*ones(8,1); 4*ones(8,1); ...
                      5*ones(8,1); 6*ones(8,1); 7*ones(8,1); 8*ones(8,1)];
             
chanListAscending = 1:64;


duration  = 0.01;
fs        = 30000;
binSize   = 1/fs;

preLength  = 36;
postLength = 36;

jscale     = 1;
njitter    = 500;
alpha      = 0.05;
for_grant  = false;
fig_use    = 102;

tWave = (1/30)*(-preLength+1:postLength);

pairsAna = [3   44;             45  20;             20  79;             34  79;
            79  91;             79  118;            91  118;            44  79;
            91  109;            73  91;             119 91;             45  4;
            45  19;             79  38;             3   45;             4   79;
            18  45;             32  79;             79  38;             3   38;
            91  101;            45  12;             82  91;             3   25;
            70  91;             20  34;             20  44;             40  79;
            79  58;             3   32;             4   25;             20  38;
            46  79;             79  87;             87  119;            98  119;
            119 122;            45  15;             17  45;             79  19;
            91  68;             79  95;             91  110;            79  45; 
            38  3];

maxSimilarityChs = zeros(size(pairsAna));
        
for loopPairs = 1:size(pairsAna)
    
    refID = pairsAna(loopPairs,1);
    tarID = pairsAna(loopPairs,2);

    refIdx = find(unitsGlobalSnapshots.neuronID == refID);
    tarIdx = find(unitsGlobalSnapshots.neuronID == tarID);

    refWave = unitsGlobalSnapshots.waveforms{refIdx};
    tarWave = unitsGlobalSnapshots.waveforms{tarIdx};

    ymin = min(min([refWave tarWave]));
    ymax = max(max([refWave tarWave]));
    
    %%
    
    [~,refPeakChan] = min(min(refWave'));
    [~,tarPeakChan] = min(min(tarWave'));
    
    refShank = ceil(refPeakChan/8);
    tarShank = ceil(tarPeakChan/8);
    
    %%
    
    refChIdxTemp = find(channelShankLayout == refShank);
    tarChIdxTemp = find(channelShankLayout == tarShank);
    
    refTemp       = min(refWave(refChIdxTemp,:),[],2);
    tarTemp       = min(tarWave(refChIdxTemp,:),[],2);
    [~,idxTemp]   = min(abs(refTemp - tarTemp));
    refMaxSimChan = refChIdxTemp(idxTemp);
    
    refTemp       = min(refWave(tarChIdxTemp,:),[],2);
    tarTemp       = min(tarWave(tarChIdxTemp,:),[],2);
    [~,idxTemp]   = min(abs(refTemp - tarTemp));
    tarMaxSimChan = tarChIdxTemp(idxTemp);
    
    maxSimilarityChs(loopPairs,1) = refMaxSimChan;
    maxSimilarityChs(loopPairs,2) = tarMaxSimChan;
    
    %%

    spike_data_fullpath = '/home/nasko/CUNY_Work_Han_Kamran_Nat/data/wake_new/wake-spikes.mat';
    dataPath = '/home/nasko/CUNY_Work_Han_Kamran_Nat/data/wake_new/';

    [data_dir, name, ~] = fileparts(spike_data_fullpath);

    load(spike_data_fullpath, 'spikes')
    load(fullfile(data_dir, 'wake-behavior.mat'), 'behavior');
    load(fullfile(data_dir, 'wake-basics.mat'),   'basics');

    temp = spikes.RoyMaze1(refID).time;
    temp = temp((temp > behavior.RoyMaze1.time(2,1)) & (temp < behavior.RoyMaze1.time(2,2)));
    refSpikeTimes = temp/1e6; % convert to ms
    
    temp = spikes.RoyMaze1(tarID).time;
    temp = temp((temp > behavior.RoyMaze1.time(2,1)) & (temp < behavior.RoyMaze1.time(2,2)));
    tarSpikeTimes = temp/1e6; % convert to ms

    %% union spike times

    for loopShank = 1:8

        refShNoiseST             = load([dataPath 'RoyMaze_shank' num2str(loopShank) '_spikeTimesNoise']);
        refShValidUnitsST        = load([dataPath 'RoyMaze_shank' num2str(loopShank) '_spikeTimesValidUnits']);
        refShValidUnitsST        = setdiff(refShValidUnitsST.validUnits,refSpikeTimes);
        refShankAllST{loopShank} = [refShNoiseST.noiseCluster; refShValidUnitsST];

    end   
    
    %% synchronous spikes
    
    figure(fig_use)
    
    [GSPExc,GSPInh,pvalE,pvalI,ccgR,tR,LSPExc,LSPInh,JBSIE,JBSII] = ...
      CCG_jitter(refSpikeTimes,tarSpikeTimes,fs,binSize,duration,'jscale',jscale, ...
                    'plot_flag', true, ...
                    'plot_output', get(fig_use, 'Number'), ...
                    'njitter', njitter, 'alpha', alpha,...
                    'for_grant', for_grant,  'plot_pointwiseBands',true);
    
    refSpikeTimes = refSpikeTimes';
    tarSpikeTimes = tarSpikeTimes';
                
    [SyncSp,SyncCCG,SyncSpBinAll] = SyncSpikes(refSpikeTimes, tarSpikeTimes, ones(length(GSPExc),1),duration);
    
    refSpikeTimesSync = [];
    tarSpikeTimesSync = [];
    
    refSpikeTimesNoSync = [];
    tarSpikeTimesNoSync = [];
    
    for k = 1:length(GSPExc)
       
        if (GSPExc(k) == 1) 
            refSpikeTimesSync = [refSpikeTimesSync; SyncSpBinAll{k}(:,1)];
            tarSpikeTimesSync = [tarSpikeTimesSync; SyncSpBinAll{k}(:,2)];
        end
            
    end
    
    
    for k = 1:length(GSPInh)
       
        if (GSPInh(k) == 1)
            refSpikeTimesNoSync = [refSpikeTimesNoSync; SyncSpBinAll{k}(:,1)];
            tarSpikeTimesNoSync = [tarSpikeTimesNoSync; SyncSpBinAll{k}(:,2)];
        end
            
    end
    
    refSpikeTimesSpon = setdiff(refSpikeTimes,[refSpikeTimesSync; refSpikeTimesNoSync]);
    tarSpikeTimesSpon = setdiff(tarSpikeTimes,[tarSpikeTimesSync; tarSpikeTimesNoSync]);
    
    refSpikeTimesSpon = refSpikeTimesSpon(randperm(length(refSpikeTimesSpon),length(refSpikeTimesSync)));
    tarSpikeTimesSpon = tarSpikeTimesSpon(randperm(length(tarSpikeTimesSpon),length(tarSpikeTimesSync)));
    
    %% synchrony evoked LFP
    
    session_name = 'RoyMaze1';
    fpass        = 300;
    
    preLength  = 3e3;
    postLength = 3e4;
    
    binSize        = 0.01*fs;
    color          = 'k';
    
    t = (-0.1+(1/fs)):1/fs:1;
    
    [chanData, onsetTime] = HiroLoadRawNSC(session_name,refPeakChan);
    
    refSpikeTimeIndex = round((refSpikeTimesSync-onsetTime/1e6)*fs) + 1;
%     refSpikeTimeIndex = refSpikeTimeIndex(randperm(length(refSpikeTimeIndex),10000));
    
    waveAvgSync = waveformAvg(double(chanData),refSpikeTimeIndex, preLength,postLength,fpass,fs,false,true);
    
    refSpikeTimeIndex = round((refSpikeTimesNoSync-onsetTime/1e6)*fs) + 1;
%     refSpikeTimeIndex = refSpikeTimeIndex(randperm(length(refSpikeTimeIndex),10000));
    
    waveAvgNoSync = waveformAvg(double(chanData), refSpikeTimeIndex, preLength,postLength,fpass,fs,false,true);
    
    refSpikeTimeIndex = round((refSpikeTimesSpon-onsetTime/1e6)*fs) + 1;
%     refSpikeTimeIndex = refSpikeTimeIndex(randperm(length(refSpikeTimeIndex),length(refSpikeTimesSync)));
    
    waveAvgSpon = waveformAvg(double(chanData), refSpikeTimeIndex, preLength,postLength,fpass,fs,false,true);
    
    tiledlayout(2,2)
    
    nexttile
    plot(t,waveAvgSync)
    hold on
%     plot(t,waveAvgNoSync)
    plot(t,waveAvgSpon)
    hold off
    title('LFP at ref peak channel')
%     legend('ref sync raw LFP','ref no sync raw LFP','ref spon raw LFP')
    legend('ref sync raw LFP','ref spon raw LFP')
    xlabel('[s]')
    ylabel('[mV]')
    xlims = get(gca,'xlim');
    xlim(xlims)
    
    [chanData, onsetTime] = HiroLoadRawNSC(session_name,tarPeakChan);
    
    tarSpikeTimeIndex = round((tarSpikeTimesSync-onsetTime/1e6)*fs) + 1;
%     tarSpikeTimeIndex = tarSpikeTimeIndex(randperm(length(tarSpikeTimeIndex),10000));
    
    waveAvgSync = waveformAvg(double(chanData),tarSpikeTimeIndex, preLength,postLength,fpass,fs,false,true);
    
    tarSpikeTimeIndex = round((tarSpikeTimesNoSync-onsetTime/1e6)*fs) + 1;
%     tarSpikeTimeIndex = tarSpikeTimeIndex(randperm(length(tarSpikeTimeIndex),10000));
    
    waveAvgNoSync = waveformAvg(double(chanData), tarSpikeTimeIndex, preLength,postLength,fpass,fs,false,true);
    
    tarSpikeTimeIndex = round((refSpikeTimesSpon-onsetTime/1e6)*fs) + 1;
%     tarSpikeTimeIndex = tarSpikeTimeIndex(randperm(length(tarSpikeTimeIndex),length(refSpikeTimesSync)));
    
    waveAvgSpon = waveformAvg(double(chanData), refSpikeTimeIndex, preLength,postLength,fpass,fs,false,true);
    
    nexttile
    plot(t,waveAvgSync)
    hold on
%     plot(t,waveAvgNoSync)
    plot(t,waveAvgSpon)
    hold off
    title('LFP at tar peak channel')
%     legend('tar sync raw LFP','tar no sync raw LFP','tar spon raw LFP')
    legend('tar sync raw LFP','tar spon raw LFP')
    xlabel('[s]')
    ylabel('[mV]')
    xlims = get(gca,'xlim');
    xlim(xlims)
    
    t = (-0.1+(1/(fs/binSize))):1/(fs/binSize):1;
    
    [refNSync,refNperTrialSync] = PSTH(refSpikeTimes,refSpikeTimesSync,postLength,preLength,binSize,fs,color);
    [refNSpon,refNperTrialSpon] = PSTH(refSpikeTimes,refSpikeTimesSpon,postLength,preLength,binSize,fs,color);
    
    nexttile
    plot(t,refNSync); hold on; plot(t,refNSpon); hold off
    legend('ref sync PSTH','ref spon PSTH')
    ylabel('counts')
    xlims = get(gca,'xlim');
    xlim(xlims)
    
    [tarNSync,tarNperTrialSync] = PSTH(tarSpikeTimes,tarSpikeTimesSync,postLength,preLength,binSize,fs,color);
    [tarNSpon,tarNperTrialSpon] = PSTH(tarSpikeTimes,tarSpikeTimesSpon,postLength,preLength,binSize,fs,color);
    
    nexttile
    plot(t,tarNSync); hold on; plot(t,tarNSpon); hold off
    legend('tar sync PSTH','tar spon PSTH')
    ylabel('counts')
    xlims = get(gca,'xlim');
    xlim(xlims)
    
    %% global snapshot

    figure
    tiledlayout(5,16)

    for i = 1:64

        nexttile

        plot(tWave,refWave(channelLayout(i),:),'LineWidth',2)
%         hold on
%         plot(tWave,tarWave(channelLayout(i),:),'LineWidth',2)
%         hold off

        ylim([ymin ymax])
        set(gca,'xtick',[])
        set(gca,'xticklabel',[])
        set(gca,'ytick',[])
        set(gca,'yticklabel',[])
        set(gca,'YDir','reverse')

    end

    for loopShank = 1:8

        [ccg,tR] = CCG([refSpikeTimes', refShankAllST{loopShank}'],...
                    [ones(size(refSpikeTimes')) 2*ones(size(refShankAllST{loopShank}'))], ...
                    'binSize', binSize, 'duration', duration, 'Fs', 1/fs, ...
                    'norm', 'counts');
        ccg=ccg(:,1,2);
        
        [ccgSynch,tR] = CCG([refSpikeTimesSync', refShankAllST{loopShank}'],...
                    [ones(size(refSpikeTimesSync')) 2*ones(size(refShankAllST{loopShank}'))], ...
                    'binSize', binSize, 'duration', duration, 'Fs', 1/fs, ...
                    'norm', 'counts');
        ccgSynch=ccgSynch(:,1,2);
        
        tR = tR*1000; % convert to ms

        nexttile
        plot(tR,ccg/sum(ccg))
        hold on
        plot(tR,ccgSynch/sum(ccgSynch))
        hold
        xlim([min(tWave) max(tWave)])
        set(gca,'xtick',[])
        set(gca,'xticklabel',[])
        set(gca,'ytick',[])
        set(gca,'yticklabel',[])

        nexttile
        plot(tR,ccg/sum(ccg))
        hold on
        plot(tR,ccgSynch/sum(ccgSynch))
        hold
        xlim([min(tWave) max(tWave)])
        set(gca,'xtick',[])
        set(gca,'xticklabel',[])
        set(gca,'ytick',[])
        set(gca,'yticklabel',[])

    end
    
    close all
    
end











