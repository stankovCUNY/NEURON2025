function runSpWaveAnaBursting(state)

    close all
    
    %% metric option:
    %
    % 'troughAmplitude'
    % 'peakToTrough'
    % 'halfWidth'
    % 'spikeWidth'

    %% state options:
    % 
    % 'allStates'
    % 'spindle'
    % 'ripple'
    % 'mobility'
    % 'theta'
    % 'gamma'
    % 'noTheta'
    % 'noGamma'
    % 'noThetaAndGamma'
    % 'noStates'

    % condition markers
    [spindleMaze,rippleMaze,mobileMaze,thetaMaze,gammaMaze] = makeBXandFreqFlags;


    %% hardcoding pairs of interest for waveform analysis

    % constants
    session_name   = 'RoyMaze1';
    % session_name   = 'KevinMaze1';
    conn_type      = {'GapPairs','ExcPairs','InhPairs'};
    jscale         = 1;
    alpha_name     = 5;
%     duration       = CCGwidth;
%     duration       = 0.002;
    duration       = 0.01;
    fpass          = 300;
    fig_use        = 102;
    njitter        = 500;
    alpha          = 0.05;
    for_grant      = false;
    filterFlag     = false;
    Nwaveforms     = 100;
    resolution_use = '-r300'; % >= 300dpi required by Nature Neuro, so that's what I'm using
    fs             = 30000;
    binSize        = 1/fs;

    preLength  = 36;
    postLength = 36;

    tWave     = (-26:27)*(1/30);
    spikeTime = (1/30)*(-preLength+1:postLength);

    if strcmp(session_name,'RoyMaze1')
%         pairsAna = [16 38]; % only monosynapse
        pairsAna = [3   44;             20  45;             20  79;             34  79;
                    79  91;             79  118;            91  118;            44  79;
                    91  109;            73  91;             119 91;             45  4;
                    45  19;             79  38;             3   45;             4   79;
                    18  45;             32  79;             79  38;             3   38;
                    91  101;            45  12;             82  91;             3   25;
                    70  91;             20  34;             20  44;             40  79;
                    79  58;             3   32;             4   25;             20  38;
                    46  79;             79  87;             87  119;            98  119;
                    119 122;            45  15;             17  45;             79  19;
                    91  68;             79  95;             91  110;            79  45];

    elseif strcmp(session_name,'KevinMaze1')

        pairsAna = [70  22;             70  54;             70  50;             3   22];

    end

    shankChanList ={1:8  , ...
                    9:16 , ...
                    17:24, ...
                    25:32, ...
                    33:40, ...
                    41:48, ...
                    49:56, ...
                    57:64};

    % Monitor specific plot settings.
    screensize = get(0,'screensize');
    % initiate figure
    hcomb = figure(102);

    res_type = 'QHD';
    pos = [70 230 2660 1860]; a_offset = [0 850 100 900]'; b_offset = [0 0 -100 -100]';
    arrayfun(@(a) set(a, 'Position', pos), hcomb(:));

    % pair data 
    datapath = '/home/nasko/CUNY_Work_Han_Kamran_Nat/data/wake_new/';

    [jit_var_out,conn_type_idx] = loadJitVarOut(session_name,conn_type,jscale,alpha_name,datapath);

    % get all pairs
    cell_pairs = {jit_var_out.cell_pair}';

    % convert to matrix
    cell_pairs_mat = [];
    for i = 1:size(cell_pairs,1)
        cell_pairs_mat = [cell_pairs_mat; cell_pairs{i}];
    end

    %% get session maze time markers
    spike_data_fullpath = "/home/nasko/CUNY_Work_Han_Kamran_Nat/data/wake_new/wake-spikes.mat";
    [data_dir, name, ~] = fileparts(spike_data_fullpath);
    load(fullfile(data_dir, 'wake-behavior.mat'), 'behavior');
    load(fullfile(data_dir, 'wake-spikes.mat'));


    %% get channels from waveform files. Not computationally savvy but it's fastest solution
    extraDrivePath = '/media/nasko/3ce6d8b6-2570-4c8f-941f-9b4b0641301c/home/nasko/CUNY_Work_Han_Kamran_Hiro_Data/';

    if strcmp(session_name,'RoyMaze1')
        Session = bz_GetSpikes_Atanas('basepath',[extraDrivePath 'Roy-maze1'],'noPrompts',true);
    elseif strcmp(session_name,'KevinMaze1')
        Session = bz_GetSpikes_Atanas('basepath',[extraDrivePath 'Kevin-maze1'],'noPrompts',true);
    end
    
    %% save variables
    burstSplit.pairs          = [];
    burstSplit.refShortPreISI = {};
    burstSplit.refLongPreISI  = {};
    burstSplit.tarShortPreISI = {};
    burstSplit.tarLongPreISI  = {};
    burstSplit.refCh          = [];
    burstSplit.tarCh          = [];
    
    %%

    waveDataSpCounts = [];
    for i = 1:size(Session.times,2)
        waveDataSpCounts(i) = size(Session.times{i},1);
    end

    for j = 24:size(pairsAna,1)
        
        tiledlayout(5,2)
        
        idx = j; % pair to run, right now we are working with pre-maze time period

        current_pair = pairsAna(idx,:);

        current_pair_indices = find((current_pair(1) == cell_pairs_mat(:,1)) & ...
                                    (current_pair(2) == cell_pairs_mat(:,2)));
        
        burstSplit.pairs = [burstSplit.pairs; current_pair];
        %%

        n1maze  = jit_var_out(current_pair_indices(2)).cell1_spike_times;
        n2maze  = jit_var_out(current_pair_indices(2)).cell2_spike_times;

        n1shank = jit_var_out(current_pair_indices(2)).cell1shank;
        n2shank = jit_var_out(current_pair_indices(2)).cell2shank;

        numSpikesCell1 = length(jit_var_out(current_pair_indices(1)).cell1_spike_times) + ...
                         length(jit_var_out(current_pair_indices(2)).cell1_spike_times) + ...
                         length(jit_var_out(current_pair_indices(3)).cell1_spike_times);

        numSpikesCell2 = length(jit_var_out(current_pair_indices(1)).cell2_spike_times) + ...
                         length(jit_var_out(current_pair_indices(2)).cell2_spike_times) + ...
                         length(jit_var_out(current_pair_indices(3)).cell2_spike_times);

        waveIdxCell1 = find(numSpikesCell1 == waveDataSpCounts);
        waveIdxCell2 = find(numSpikesCell2 == waveDataSpCounts);

        chCell1 = Session.maxWaveformCh(waveIdxCell1);
        chCell2 = Session.maxWaveformCh(waveIdxCell2);

        cell1type = jit_var_out(current_pair_indices(2)).cell1type;
        cell2type = jit_var_out(current_pair_indices(2)).cell2type;

        cellID1 = current_pair(1);
        cellID2 = current_pair(2);

        % plot CCG 
        res1 = jit_var_out(current_pair_indices(2)).cell1_spike_times;
        res2 = jit_var_out(current_pair_indices(2)).cell2_spike_times;

        if strcmp(state,'spindle')

            res1 = filterSpikes(res1,spindleMaze);
            res2 = filterSpikes(res2,spindleMaze);

        elseif strcmp(state,'ripple')

            res1 = filterSpikes(res1,rippleMaze);
            res2 = filterSpikes(res2,rippleMaze);

        elseif strcmp(state,'mobility')

            res1 = filterSpikes(res1,mobileMaze);
            res2 = filterSpikes(res2,mobileMaze);

        elseif strcmp(state,'theta')

            res1 = filterSpikes(res1,thetaMaze);
            res2 = filterSpikes(res2,thetaMaze);

        elseif strcmp(state,'gamma')

            res1 = filterSpikes(res1,gammaMaze);
            res2 = filterSpikes(res2,gammaMaze);
        
        elseif strcmp(state,'noRipple')
                        
            res1Ripple = filterSpikes(res1,rippleMaze);
            res1       = setdiff(res1,res1Ripple);

            res2Ripple = filterSpikes(res2,rippleMaze);
            res2       = setdiff(res2,res2Ripple);
            
        elseif strcmp(state,'noTheta')

            res1Theta = filterSpikes(res1,thetaMaze);
            res1      = setdiff(res1,res1Theta);

            res2Theta = filterSpikes(res2,thetaMaze);
            res2      = setdiff(res2,res2Theta);

        elseif strcmp(state,'noGamma')

            res1Gamma = filterSpikes(res1,gammaMaze);
            res1      = setdiff(res1,res1Gamma);

            res2Gamma = filterSpikes(res2,gammaMaze);
            res2      = setdiff(res2,res2Gamma);
        
        elseif strcmp(state,'noThetaAndRipple')
            
            res1Theta    = filterSpikes(res1,thetaMaze);
            res1Ripple   = filterSpikes(res1,rippleMaze);

            res1ThetaAndRipple = [res1Theta; res1Ripple];
            res1ThetaAndRipple = sort(unique(res1ThetaAndRipple));
            res1               = setdiff(res1,res1ThetaAndRipple);

            res2Theta    = filterSpikes(res2,thetaMaze);
            res2Ripple   = filterSpikes(res2,rippleMaze);

            res2ThetaAndRipple = [res2Theta; res2Ripple];
            res2ThetaAndRipple = sort(unique(res2ThetaAndRipple));
            res2               = setdiff(res2,res2ThetaAndRipple);
            
        elseif strcmp(state,'noThetaAndGamma')

            res1Theta    = filterSpikes(res1,thetaMaze);
            res1Gamma    = filterSpikes(res1,gammaMaze);

            res1ThetaAndRipple = [res1Theta; res1Gamma];
            res1ThetaAndRipple = sort(unique(res1ThetaAndRipple));
            res1               = setdiff(res1,res1ThetaAndRipple);

            res2Theta    = filterSpikes(res2,thetaMaze);
            res2Gamma    = filterSpikes(res2,gammaMaze);

            res2ThetaAndRipple = [res2Theta; res2Gamma];
            res2ThetaAndRipple = sort(unique(res2ThetaAndRipple));
            res2               = setdiff(res2,res2ThetaAndRipple);

        elseif strcmp(state,'noStates')

            res1Spindle  = filterSpikes(res1,spindleMaze);
            res1Ripple   = filterSpikes(res1,rippleMaze);
            res1Mobility = filterSpikes(res1,mobileMaze);
            res1Theta    = filterSpikes(res1,thetaMaze);
            res1Gamma    = filterSpikes(res1,gammaMaze);

            res1State = [res1Spindle; res1Ripple; res1Mobility; res1Theta; res1Gamma];
            res1State = sort(unique(res1State));
            res1      = setdiff(res1,res1State);

            res2Spindle  = filterSpikes(res2,spindleMaze);
            res2Ripple   = filterSpikes(res2,rippleMaze);
            res2Mobility = filterSpikes(res2,mobileMaze);
            res2Theta    = filterSpikes(res2,thetaMaze);
            res2Gamma    = filterSpikes(res2,gammaMaze);

            res2State = [res2Spindle; res2Ripple; res2Mobility; res2Theta; res2Gamma];
            res2State = sort(unique(res2State));
            res2      = setdiff(res2,res2State);

        end
        
        [chanDataAtCell1, onsetTime] = HiroLoad300hz(session_name,chCell1);

        pairStr = [num2str(cellID1) cell1type ' (ch: ' num2str(chCell1) ', sh: ' num2str(n1shank) ')' ' v ' ...
                   num2str(cellID2) cell2type ' (ch: ' num2str(chCell2) ', sh: ' num2str(n2shank) ')'];

        saveStr = [session_name ' (burst split - condition: ' state ') - '  num2str(cellID1) cell1type ' (ch: ' num2str(chCell1) ', sh: ' num2str(n1shank) ')' ' v ' ...
                                                                             num2str(cellID2) cell2type ' (ch: ' num2str(chCell2) ', sh: ' num2str(n2shank) ')'];

        % get waveforms
        spikeTimeIndxCell1 = round((res1-onsetTime/1e6)*fs) + 1; 

        [spikeAvgMaxChanCell1, waveformsMaxCell1] = waveformAvg(chanDataAtCell1,spikeTimeIndxCell1,preLength,postLength,fpass,fs,filterFlag);    

%         title(['eSpike: ' num2str(cellID1) cell1type ' (ch: ' num2str(chCell1) ', sh: ' num2str(n1shank) ') condition: ' state])
        
        ISI    = diff(res1);
        ISIbin = ISI < 0.004;
        ISIbin = [false; ISIbin]; % correct for lost sample
             
        resShortPreISI  = res1(ISIbin);
        resLongPreISI   = res1(~ISIbin);
%         resLongPreISI   = resLongPreISI(randperm(length(resLongPreISI),length(resShortPreISI)));
        
        burstSplit.refShortPreISI{j} = resShortPreISI;
        burstSplit.refLongPreISI{j}  = resLongPreISI;
        burstSplit.refCh          = [burstSplit.refCh chCell1];
        burstSplit.tarCh          = [burstSplit.tarCh chCell2];
    
        waveShortPreISI = waveformsMaxCell1(:,ISIbin);
        waveLongPreISI  = waveformsMaxCell1(:,~ISIbin);
        
        if  size(waveShortPreISI,2) < Nwaveforms
            waveShortPreISIRandIdx = 1:size(waveShortPreISI,2);
        else
            waveShortPreISIRandIdx = randperm(size(waveShortPreISI,2),Nwaveforms);
        end
        waveLongPreISIRandIdx  = randperm(size(waveLongPreISI,2),Nwaveforms);

        hcomb = figure(102);
        [GSPExc,GSPInh,pvalE,pvalI,ccgR,tR,LSPExc,LSPInh,JBSIE,JBSII] = ...
          CCG_jitter(res1,res2,fs,binSize,duration,'jscale',jscale, ...
                        'plot_flag', false, ...
                        'plot_output', get(fig_use, 'Number'), ...
                        'njitter', njitter, 'alpha', alpha,...
                        'for_grant', for_grant,  'plot_pointwiseBands',false);

        hcomb = figure(102);
        [GSPExc,GSPInh,pvalE,pvalI,ccgRShortPreISI,tR,LSPExc,LSPInh,JBSIE,JBSII] = ...
          CCG_jitter(resShortPreISI,res2,fs,binSize,duration,'jscale',jscale, ...
                        'plot_flag', false, ...
                        'plot_output', get(fig_use, 'Number'), ...
                        'njitter', njitter, 'alpha', alpha,...
                        'for_grant', for_grant,  'plot_pointwiseBands',false);

        hcomb = figure(102);
        [GSPExc,GSPInh,pvalE,pvalI,ccgRLongPreISI,tR,LSPExc,LSPInh,JBSIE,JBSII] = ...
          CCG_jitter(resLongPreISI,res2,fs,binSize,duration,'jscale',jscale, ...
                        'plot_flag', false, ...
                        'plot_output', get(fig_use, 'Number'), ...
                        'njitter', njitter, 'alpha', alpha,...
                        'for_grant', for_grant,  'plot_pointwiseBands',false);

        nexttile(1)
        patchline(spikeTime,waveShortPreISI(:,waveShortPreISIRandIdx(1)),'EdgeColor','#D95319','LineWidth',0.25,'EdgeAlpha',1);
        hold on
        for k = 2:length(waveShortPreISIRandIdx)
            patchline(spikeTime,waveShortPreISI(:,waveShortPreISIRandIdx(k)),'EdgeColor','#D95319','LineWidth',0.25,'EdgeAlpha',1);
        end
        hold off
        set(gca,'YDir','reverse')
        title(['short ISIs (<4ms, n = ' num2str(size(waveShortPreISI,2)) ')'])
        
        nexttile(3)
        patchline(spikeTime,waveLongPreISI(:,waveLongPreISIRandIdx(1)),'EdgeColor','#EDB120','LineWidth',0.25,'EdgeAlpha',1);
        hold on
        for k = 2:Nwaveforms
            patchline(spikeTime,waveLongPreISI(:,waveLongPreISIRandIdx(k)),'EdgeColor','#EDB120','LineWidth',0.25,'EdgeAlpha',1);
        end
        hold off
        set(gca,'YDir','reverse')
        title(['long ISIs (>4ms, n = ' num2str(size(waveLongPreISI,2)) ')'])
        
        nexttile(5)
        plot(spikeTime,spikeAvgMaxChanCell1,'LineWidth',2)
        hold on
        plot(spikeTime,mean(waveShortPreISI,2),'LineWidth',2);
        plot(spikeTime,mean(waveLongPreISI,2),'LineWidth',2);
        hold off
        set(gca,'YDir','reverse')
        legend(['all (n = ' num2str(size(waveformsMaxCell1,2)) ')'], ...
               ['short ISIs (<4ms, n = ' num2str(size(waveShortPreISI,2)) ')'], ...
               ['long ISIs (>4ms, n = '  num2str(size(waveLongPreISI,2)) ')'])

        nexttile(7)
        plot(tR*1000,ccgR(:,1,2)/sum(ccgR(:,1,2)),'LineWidth',2)
        hold on
        plot(tR*1000,ccgRShortPreISI(:,1,2)/sum(ccgRShortPreISI(:,1,2)),'LineWidth',2)
        plot(tR*1000,ccgRLongPreISI(:,1,2)/sum(ccgRLongPreISI(:,1,2)),'LineWidth',2)
        hold off
        title(['normalized CCG: ' pairStr])
        legend('all','short ISIs (<4ms)','long ISIs (>4ms)')
        
        nexttile(9)
        plot(tR*1000,ccgR(:,1,1)/sum(ccgR(:,1,1)),'LineWidth',2)
        hold on
        plot(tR*1000,ccgRShortPreISI(:,1,1)/sum(ccgRShortPreISI(:,1,1)),'LineWidth',2)
        plot(tR*1000,ccgRLongPreISI(:,1,1)/sum(ccgRLongPreISI(:,1,1)),'LineWidth',2)
        hold off
        title(['normalized ACG: ' num2str(cellID1) cell1type ' (ch: ' num2str(chCell1) ', sh: ' num2str(n1shank) ')'])
        legend('all','short ISIs (<4ms)','long ISIs (>4ms)')

        %%

        n2maze  = jit_var_out(current_pair_indices(2)).cell1_spike_times;
        n1maze  = jit_var_out(current_pair_indices(2)).cell2_spike_times;

        n2shank = jit_var_out(current_pair_indices(2)).cell1shank;
        n1shank = jit_var_out(current_pair_indices(2)).cell2shank;

        numSpikesCell2 = length(jit_var_out(current_pair_indices(1)).cell1_spike_times) + ...
                         length(jit_var_out(current_pair_indices(2)).cell1_spike_times) + ...
                         length(jit_var_out(current_pair_indices(3)).cell1_spike_times);

        numSpikesCell1 = length(jit_var_out(current_pair_indices(1)).cell2_spike_times) + ...
                         length(jit_var_out(current_pair_indices(2)).cell2_spike_times) + ...
                         length(jit_var_out(current_pair_indices(3)).cell2_spike_times);

        waveIdxCell2 = find(numSpikesCell1 == waveDataSpCounts);
        waveIdxCell1 = find(numSpikesCell2 == waveDataSpCounts);

        chCell2 = Session.maxWaveformCh(waveIdxCell1);
        chCell1 = Session.maxWaveformCh(waveIdxCell2);

        cell2type = jit_var_out(current_pair_indices(2)).cell1type;
        cell1type = jit_var_out(current_pair_indices(2)).cell2type;

        cellID2 = current_pair(1);
        cellID1 = current_pair(2);

        % plot CCG 
        res2 = jit_var_out(current_pair_indices(2)).cell1_spike_times;
        res1 = jit_var_out(current_pair_indices(2)).cell2_spike_times;

        if strcmp(state,'spindle')

            res1 = filterSpikes(res1,spindleMaze);
            res2 = filterSpikes(res2,spindleMaze);

        elseif strcmp(state,'ripple')

            res1 = filterSpikes(res1,rippleMaze);
            res2 = filterSpikes(res2,rippleMaze);

        elseif strcmp(state,'mobility')

            res1 = filterSpikes(res1,mobileMaze);
            res2 = filterSpikes(res2,mobileMaze);

        elseif strcmp(state,'theta')

            res1 = filterSpikes(res1,thetaMaze);
            res2 = filterSpikes(res2,thetaMaze);

        elseif strcmp(state,'gamma')

            res1 = filterSpikes(res1,gammaMaze);
            res2 = filterSpikes(res2,gammaMaze);

        elseif strcmp(state,'noTheta')

            res1Theta = filterSpikes(res1,thetaMaze);
            res1      = setdiff(res1,res1Theta);

            res2Theta = filterSpikes(res2,thetaMaze);
            res2      = setdiff(res2,res2Theta);

        elseif strcmp(state,'noGamma')

            res1Gamma = filterSpikes(res1,gammaMaze);
            res1      = setdiff(res1,res1Gamma);

            res2Gamma = filterSpikes(res2,gammaMaze);
            res2      = setdiff(res2,res2Gamma);

        elseif strcmp(state,'noThetaAndGamma')

            res1Theta    = filterSpikes(res1,thetaMaze);
            res1Gamma    = filterSpikes(res1,gammaMaze);

            res1ThetaAndRipple = [res1Theta; res1Gamma];
            res1ThetaAndRipple = sort(unique(res1ThetaAndRipple));
            res1      = setdiff(res1,res1ThetaAndRipple);

            res2Theta    = filterSpikes(res2,thetaMaze);
            res2Gamma    = filterSpikes(res2,gammaMaze);

            res2ThetaAndRipple = [res2Theta; res2Gamma];
            res2ThetaAndRipple = sort(unique(res2ThetaAndRipple));
            res2      = setdiff(res2,res2ThetaAndRipple);

        elseif strcmp(state,'noStates')

            res1Spindle  = filterSpikes(res1,spindleMaze);
            res1Ripple   = filterSpikes(res1,rippleMaze);
            res1Mobility = filterSpikes(res1,mobileMaze);
            res1Theta    = filterSpikes(res1,thetaMaze);
            res1Gamma    = filterSpikes(res1,gammaMaze);

            res1State = [res1Spindle; res1Ripple; res1Mobility; res1Theta; res1Gamma];
            res1State = sort(unique(res1State));
            res1      = setdiff(res1,res1State);

            res2Spindle  = filterSpikes(res2,spindleMaze);
            res2Ripple   = filterSpikes(res2,rippleMaze);
            res2Mobility = filterSpikes(res2,mobileMaze);
            res2Theta    = filterSpikes(res2,thetaMaze);
            res2Gamma    = filterSpikes(res2,gammaMaze);

            res2State = [res2Spindle; res2Ripple; res2Mobility; res2Theta; res2Gamma];
            res2State = sort(unique(res2State));
            res2      = setdiff(res2,res2State);

        end
        
        [chanDataAtCell1, onsetTime] = HiroLoad300hz(session_name,chCell1);

        pairStr = [num2str(cellID1) cell1type ' (ch: ' num2str(chCell1) ', sh: ' num2str(n1shank) ')' ' v ' ...
                   num2str(cellID2) cell2type ' (ch: ' num2str(chCell2) ', sh: ' num2str(n2shank) ')'];

        res1sync = [];
        res2sync = [];

        % get waveforms
        spikeTimeIndxCell1 = round((res1-onsetTime/1e6)*fs) + 1; 
        [spikeAvgMaxChanCell1, waveformsMaxCell1] = waveformAvg(chanDataAtCell1,spikeTimeIndxCell1,preLength,postLength,fpass,fs,filterFlag);    

%         title(['eSpike: ' num2str(cellID1) cell1type ' (ch: ' num2str(chCell1) ', sh: ' num2str(n1shank) ') condition: ' state])

        ISI    = diff(res1);
        ISIbin = ISI < 0.004;
        ISIbin = [false; ISIbin]; % correct for lost sample
             
        resShortPreISI  = res1(ISIbin);
        resLongPreISI   = res1(~ISIbin);
%         resLongPreISI   = resLongPreISI(randperm(length(resLongPreISI),length(resShortPreISI)));
        
        burstSplit.refShortPreISI{j} = resShortPreISI;
        burstSplit.refLongPreISI{j}  = resLongPreISI;
        burstSplit.refCh          = [burstSplit.refCh chCell1];
        burstSplit.tarCh          = [burstSplit.tarCh chCell2];
        
        waveShortPreISI = waveformsMaxCell1(:,ISIbin);
        waveLongPreISI  = waveformsMaxCell1(:,~ISIbin);
        
        if size(waveShortPreISI,2) < Nwaveforms
            waveShortPreISIRandIdx = 1:size(waveShortPreISI,2);;
        else
            waveShortPreISIRandIdx = randperm(size(waveShortPreISI,2),Nwaveforms);
        end
        waveLongPreISIRandIdx  = randperm(size(waveLongPreISI,2),Nwaveforms);
        
        hcomb = figure(102);
        [GSPExc,GSPInh,pvalE,pvalI,ccgR,tR,LSPExc,LSPInh,JBSIE,JBSII] = ...
          CCG_jitter(res1,res2,fs,binSize,duration,'jscale',jscale, ...
                        'plot_flag', false, ...
                        'plot_output', get(fig_use, 'Number'), ...
                        'njitter', njitter, 'alpha', alpha,...
                        'for_grant', for_grant,  'plot_pointwiseBands',false);

        hcomb = figure(102);
        [GSPExc,GSPInh,pvalE,pvalI,ccgRShortPreISI,tR,LSPExc,LSPInh,JBSIE,JBSII] = ...
          CCG_jitter(resShortPreISI,res2,fs,binSize,duration,'jscale',jscale, ...
                        'plot_flag', false, ...
                        'plot_output', get(fig_use, 'Number'), ...
                        'njitter', njitter, 'alpha', alpha,...
                        'for_grant', for_grant,  'plot_pointwiseBands',false);

        hcomb = figure(102);
        [GSPExc,GSPInh,pvalE,pvalI,ccgRLongPreISI,tR,LSPExc,LSPInh,JBSIE,JBSII] = ...
          CCG_jitter(resLongPreISI,res2,fs,binSize,duration,'jscale',jscale, ...
                        'plot_flag', false, ...
                        'plot_output', get(fig_use, 'Number'), ...
                        'njitter', njitter, 'alpha', alpha,...
                        'for_grant', for_grant,  'plot_pointwiseBands',false);
        
        nexttile(2)
        patchline(spikeTime,waveShortPreISI(:,waveShortPreISIRandIdx(1)),'EdgeColor','#D95319','LineWidth',0.25,'EdgeAlpha',1);
        hold on
        for k = 2:length(waveShortPreISIRandIdx)
            patchline(spikeTime,waveShortPreISI(:,waveShortPreISIRandIdx(k)),'EdgeColor','#D95319','LineWidth',0.25,'EdgeAlpha',1);
        end
        hold off
        set(gca,'YDir','reverse')
        title(['short ISIs (<4ms, n = ' num2str(size(waveShortPreISI,2)) ')'])
        
        nexttile(4)
        patchline(spikeTime,waveLongPreISI(:,waveLongPreISIRandIdx(1)),'EdgeColor','#EDB120','LineWidth',0.25,'EdgeAlpha',1);
        hold on
        for k = 2:Nwaveforms
            patchline(spikeTime,waveLongPreISI(:,waveLongPreISIRandIdx(k)),'EdgeColor','#EDB120','LineWidth',0.25,'EdgeAlpha',1);
        end
        hold off
        set(gca,'YDir','reverse')
        title(['long ISIs (>4ms, n = ' num2str(size(waveLongPreISI,2)) ')'])
                    
        nexttile(6)
        plot(spikeTime,spikeAvgMaxChanCell1,'LineWidth',2)
        hold on
        plot(spikeTime,mean(waveShortPreISI,2),'LineWidth',2);
        plot(spikeTime,mean(waveLongPreISI,2),'LineWidth',2);
        hold off
        set(gca,'YDir','reverse')
        legend(['all (n = ' num2str(size(waveformsMaxCell1,2)) ')'], ...
               ['short ISIs (<4ms, n = ' num2str(size(waveShortPreISI,2)) ')'], ...
               ['long ISIs (>4ms, n = '  num2str(size(waveLongPreISI,2)) ')'])

        nexttile(8)
        plot(tR*1000,ccgR(:,1,2)/sum(ccgR(:,1,2)),'LineWidth',2)
        hold on
        plot(tR*1000,ccgRShortPreISI(:,1,2)/sum(ccgRShortPreISI(:,1,2)),'LineWidth',2)
        plot(tR*1000,ccgRLongPreISI(:,1,2)/sum(ccgRLongPreISI(:,1,2)),'LineWidth',2)
        hold off
        title(['normalized CCG: ' pairStr])
        legend('all','short ISIs (<4ms)','long ISIs (>4ms)')

        nexttile(10)
        plot(tR*1000,ccgR(:,1,1)/sum(ccgR(:,1,1)),'LineWidth',2)
        hold on
        plot(tR*1000,ccgRShortPreISI(:,1,1)/sum(ccgRShortPreISI(:,1,1)),'LineWidth',2)
        plot(tR*1000,ccgRLongPreISI(:,1,1)/sum(ccgRLongPreISI(:,1,1)),'LineWidth',2)
        hold off
        title(['normalized ACG: ' num2str(cellID1) cell1type ' (ch: ' num2str(chCell1) ', sh: ' num2str(n1shank) ')'])
        legend('all','short ISIs (<4ms)','long ISIs (>4ms)')
        
        %%

        save_file = fullfile(data_dir, saveStr);
        print(fig_use, save_file,'-dpng',resolution_use);

        % reset plot
        close 102

        hcomb = figure(102);
        arrayfun(@(a) set(a, 'Position', pos), hcomb(:));        
        
    end
    
    close 102

    %% save burst split subsets
    save([datapath 'burstSplitSets_' state],'burstSplit')
    
end

                         
                         
                         