function zeroLagArtifact(anaVar)
   
    pairType      = 'exquisite';
    filteredPairs = BapunScreenJitter('U',pairType);
        
    neurons    = load('/media/nasko/WD_BLACK3/RatU_temp2/RatU_Day2NSD_2021-07-24_08-16-38.neurons.mat');
    paradigm   = load('/media/nasko/WD_BLACK3/RatU_temp2/RatU_Day2NSD_2021-07-24_08-16-38.paradigm.mat');
    probegroup = load('/media/nasko/WD_BLACK3/RatU_temp2/RatU_Day2NSD_2021-07-24_08-16-38.probegroup.mat');

    % hard coding from probe file
    channelSet = {1:16, ...
                  17:32, ... 
                  33:48, ...
                  49:64, ...
                  65:80, ...
                  81:96, ...
                  97:112, ...
                  113:128, ...
                  129:144, ...
                  145:160, ...
                  161:176, ...
                  177:192};        
              
    matWaveIdx = find(cell2mat(probegroup.connected));
    
    %%

    jscale         = 5;
    alpha_name     = 5;
    duration       = 0.007;
    fs             = 30000;
    fpass          = 300;
    binSize        = 1/fs;
    fig_use        = 102;
    njitter        = 500;
    alpha          = 0.05;
    for_grant      = false;
    plotFlag       = true;
    resolution_use = '-r300'; % >= 300dpi required by Nature Neuro, so that's what I'm using.
    Nwaveforms     = 100;
    filterFlag     = true;
    
    tWave          = (-26:27)*(1/30);
    
    preLength  = 36;
    postLength = 36;
    
    samples = double(paradigm.stop(4)*neurons.sampling_rate);
    
    figpath  = '/media/nasko/3ce6d8b6-2570-4c8f-941f-9b4b0641301c/home/nasko/CUNY_Work_Han_Kamran_Utku_dataV2/figures';
    UtkuPath = '/media/nasko/3ce6d8b6-2570-4c8f-941f-9b4b0641301c/home/nasko/CUNY_Work_Han_Kamran_Utku_dataV2/data/';

    % Monitor specific plot settings.
    screensize = get(0,'screensize');
    % initiate figure
    hcomb = figure(102);

    res_type = 'QHD';
    pos = [70 230 1920 1080]; a_offset = [0 850 100 900]'; b_offset = [0 0 -100 -100]';
    arrayfun(@(a) set(a, 'Position', pos), hcomb(:));

    for i = 38:size(filteredPairs,1)
        
%         load('/media/nasko/3ce6d8b6-2570-4c8f-941f-9b4b0641301c/home/nasko/filterResRatU');
%         filterRes = double(filterRes);
%         
%         filterResBin = logical(zeros(samples,1));
%         filterResBin(filterRes) = 1;
        
        res1    = neurons.spiketrains{1,filteredPairs(i,1)};
        res2    = neurons.spiketrains{1,filteredPairs(i,2)};
        
%         res1    = filterPseudoEMG(res1);
%         res2    = filterPseudoEMG(res2);
        
%         res1idx = round(res1*fs);
%         res2idx = round(res2*fs);
%         
%         res1idx(res1idx > samples) = [];
%         res2idx(res2idx > samples) = [];
%         
%         res1bin = zeros(samples,1);
%         res2bin = zeros(samples,1); 
%         
%         res1bin(res1idx) = 1;
%         res2bin(res2idx) = 1;
%       
%         res1bin = res1bin & filterResBin;
%         res2bin = res2bin & filterResBin;
%         
%         res1 = find(res1bin)/fs;
%         res2 = find(res2bin)/fs;
        
%         res1    = setdiff(res1*fs,filterRes)/fs;
%         res2    = setdiff(res2*fs,filterRes)/fs;
        
        cellID1 = filteredPairs(i,1);
        cellID2 = filteredPairs(i,2);

        ch1     = double(neurons.peak_channels(filteredPairs(i,1))) + 1;
        ch2     = double(neurons.peak_channels(filteredPairs(i,2))) + 1;

        shank1  = double(neurons.shank_ids(filteredPairs(i,1))) + 1;
        shank2  = double(neurons.shank_ids(filteredPairs(i,2))) + 1;
        
        cell1type = neurons.neuron_type(filteredPairs(i,1),:);
        cell2type = neurons.neuron_type(filteredPairs(i,2),:);
        
        loopChan1 = channelSet{shank1};
        loopChan2 = channelSet{shank2};
        
        % reject MUAs at this stage
        if strcmp(cell1type,'mua  ') || strcmp(cell2type,'mua  ')
            continue
        end
        
        if strcmp(cell1type,'pyr  ')
            cell1type = 'p';
        elseif strcmp(cell1type,'inter')
            cell1type = 'i';
        end
        
        if strcmp(cell2type,'pyr  ')
            cell2type = 'p';
        elseif strcmp(cell2type,'inter') 
            cell2type = 'i';
        end
        
        %% 
        
%         artifactSpikeTimes = zeroLagArtifactFix;
        
%         res1 = setdiff(res1,artifactSpikeTimes);
%         res2 = setdiff(res2,artifactSpikeTimes);
        
        %% test plots 
        [GSPExc,GSPInh,pvalE,pvalI,ccgR,tR,LSPExc,LSPInh,JBSIE,JBSII] = ...
                  CCG_jitter(res1,res2,fs,binSize,duration,'jscale',jscale, 'plot_flag', plotFlag, ...
                            'plot_output', get(fig_use, 'Number'), ...
                            'njitter', njitter, 'alpha', alpha,...
                            'for_grant', for_grant);

        
        zeroLagBinOfInterest = zeros(211,1);
        zeroLagBinOfInterest(106) = 1;
        [SyncSp,SyncCCG,SyncSpBinAll,SyncSpIdxBin] = SyncSpikes(res1',res2',zeroLagBinOfInterest,duration);
        
%         plot(tR,ccgR())
        
%         res1nosync = setdiff(res1,SyncSp(:,1));
        if strcmp(anaVar,'specificShankLevel')
            
            specificShank = 1;
            
            loopChan1 = channelSet{specificShank};
            
            for j = 1:length(loopChan1)
                j
                load(['/media/nasko/WD_BLACK3/RatU_temp/highpassFiltered/neuronSnippetsPerShank/' 'rat' 'U' 'neuron' num2str(cellID1) cell1type 'shank' num2str(specificShank) 'ch' num2str(loopChan1(j)) '.mat' ])
                spikeTimesSet1{j} = data.spikeSnippets;
                clear data
            end
            
            for j = 1:length(loopChan2)
                j
                load(['/media/nasko/WD_BLACK3/RatU_temp/highpassFiltered/neuronSnippetsPerShank/' 'rat' 'U' 'neuron' num2str(cellID2) cell2type 'shank' num2str(specificShank) 'ch' num2str(loopChan2(j)) '.mat' ])
                spikeTimesSet2{j} = data.spikeSnippets;
                clear data
            end
            
            artifactIdx1 = SyncSpIdxBin{106}(:,1); artifactIdx1(artifactIdx1 > size(spikeTimesSet1{1},2)) = [];
            artifactIdx2 = SyncSpIdxBin{106}(:,2); artifactIdx2(artifactIdx2 > size(spikeTimesSet2{1},2)) = [];
            
            nonArtifactIdx1 = 1:length(res1); nonArtifactIdx1(nonArtifactIdx1 > size(spikeTimesSet1{1},2)) = []; nonArtifactIdx1(artifactIdx1) = [];
            nonArtifactIdx2 = 1:length(res2); nonArtifactIdx2(nonArtifactIdx2 > size(spikeTimesSet2{1},2)) = []; nonArtifactIdx2(artifactIdx2) = []; 
            
            for j = 1:length(loopChan1)
                nonArtifactAvg1(:,j) = mean(spikeTimesSet1{j}(:,nonArtifactIdx1),2);
                artifactAvg1(:,j)    = mean(spikeTimesSet1{j}(:,artifactIdx1),2);
            end
            
            for j = 1:length(loopChan2)
                nonArtifactAvg2(:,j) = mean(spikeTimesSet2{j}(:,nonArtifactIdx2),2);
                artifactAvg2(:,j)    = mean(spikeTimesSet2{j}(:,artifactIdx2),2);
            end
            
            figure()
            tiledlayout(3,length(loopChan1)/2)
            for j = 1:length(loopChan1)
                nexttile
                plot(nonArtifactAvg1(:,j))
                hold on
                plot(artifactAvg1(:,j))
                hold off
                title(['shank: ' num2str(shank1) ', chan: ' num2str(loopChan1(j))])
                legend non0lag 0lag 

                nexttile
                plot(nonArtifactAvg2(:,j))
                hold on
                plot(artifactAvg2(:,j))
                hold off
                title(['shank: ' num2str(shank2) ', chan: ' num2str(loopChan2(j))])
                legend non0lag 0lag 
            end
            
        elseif strcmp(anaVar,'shankLevel')

            for j = 1: length(loopChan1)
                j
                load(['/media/nasko/WD_BLACK3/RatU_temp/highpassFiltered/neuronSnippetsPerShank/' 'rat' 'U' 'neuron' num2str(cellID1) cell1type 'shank' num2str(shank1) 'ch' num2str(loopChan1(j)) '.mat' ])
                spikeTimesSet1{j} = data.spikeSnippets;
                clear data
            end
            
            for j = 1:length(loopChan2)
                j
                load(['/media/nasko/WD_BLACK3/RatU_temp/highpassFiltered/neuronSnippetsPerShank/' 'rat' 'U' 'neuron' num2str(cellID2) cell2type 'shank' num2str(shank2) 'ch' num2str(loopChan2(j)) '.mat' ])
                spikeTimesSet2{j} = data.spikeSnippets;
                clear data
            end
            
            artifactIdx1 = SyncSpIdxBin{106}(:,1); artifactIdx1(artifactIdx1 > size(spikeTimesSet1{1},2)) = [];
            artifactIdx2 = SyncSpIdxBin{106}(:,2); artifactIdx2(artifactIdx2 > size(spikeTimesSet2{1},2)) = [];
            
            nonArtifactIdx1 = 1:length(res1); nonArtifactIdx1(nonArtifactIdx1 > size(spikeTimesSet1{1},2)) = []; nonArtifactIdx1(artifactIdx1) = [];
            nonArtifactIdx2 = 1:length(res2); nonArtifactIdx2(nonArtifactIdx2 > size(spikeTimesSet2{1},2)) = []; nonArtifactIdx2(artifactIdx2) = []; 
            
            for j = 1:length(loopChan1)
                nonArtifactAvg1(:,j) = mean(spikeTimesSet1{j}(:,nonArtifactIdx1),2);
                artifactAvg1(:,j)    = mean(spikeTimesSet1{j}(:,artifactIdx1),2);
            end
            
            for j = 1:length(loopChan2)
                nonArtifactAvg2(:,j) = mean(spikeTimesSet2{j}(:,nonArtifactIdx2),2);
                artifactAvg2(:,j)    = mean(spikeTimesSet2{j}(:,artifactIdx2),2);
            end
            
            figure()
            tiledlayout(length(loopChan1),1)
            for j = 1:length(loopChan1)
                nexttile
                plot(nonArtifactAvg1(:,j))
                hold on
                plot(artifactAvg1(:,j))
                hold off
                title(['shank: ' num2str(shank1) ', chan: ' num2str(loopChan1(j))])
                legend non0lag 0lag 
                ylim([-0.5 0.4])
            end
            
            figure()
            tiledlayout(length(loopChan2),1)
            for j = 1:length(loopChan2)
                nexttile
                plot(nonArtifactAvg2(:,j))
                hold on
                plot(artifactAvg2(:,j))
                hold off
                title(['shank: ' num2str(shank2) ', chan: ' num2str(loopChan2(j))])
                legend non0lag 0lag 
                ylim([-1 0.5])
            end
        
        elseif strcmp(anaVar,'arrayLevel')
            
            channelsLoopAll = 1:192;
            channelShanks = [  ones(16,1); ...
                             2*ones(16,1); ...
                             3*ones(16,1); ...
                             4*ones(16,1); ...
                             5*ones(16,1); ...
                             6*ones(16,1);  ...
                             7*ones(16,1); ...
                             8*ones(16,1); ...
                             9*ones(16,1); ...
                             10*ones(16,1); ...
                             11*ones(16,1); ...
                             12*ones(16,1)];
                         
            for j = 1:length(channelsLoopAll)
                
                tic
                
                j
                
                load(['/media/nasko/WD_BLACK3/RatU_temp/highpassFiltered/neuronSnippetsPerShank/' 'rat' 'U' 'neuron' num2str(cellID1) cell1type 'shank' num2str(channelShanks(j)) 'ch' num2str(channelsLoopAll(j)) '.mat' ])
                
                artifactIdx1 = SyncSpIdxBin{106}(:,1); artifactIdx1(artifactIdx1 > size(data.spikeSnippets,2)) = [];
                nonArtifactIdx1 = 1:length(res1); nonArtifactIdx1(artifactIdx1) = []; nonArtifactIdx1(nonArtifactIdx1 > size(data.spikeSnippets,2)) = [];

                nonArtifactAvg1(:,j) = mean(data.spikeSnippets(:,nonArtifactIdx1),2);
                artifactAvg1(:,j)    = mean(data.spikeSnippets(:,artifactIdx1),2);
                
                clear data
                
                load(['/media/nasko/WD_BLACK3/RatU_temp/highpassFiltered/neuronSnippetsPerShank/' 'rat' 'U' 'neuron' num2str(cellID2) cell2type 'shank' num2str(channelShanks(j)) 'ch' num2str(channelsLoopAll(j)) '.mat' ])
                
                artifactIdx2 = SyncSpIdxBin{106}(:,2); artifactIdx2(artifactIdx2 > size(data.spikeSnippets,2)) = [];
                nonArtifactIdx2 = 1:length(res2); nonArtifactIdx2(artifactIdx2) = []; nonArtifactIdx2(nonArtifactIdx2 > size(data.spikeSnippets,2)) = [];

                nonArtifactAvg2(:,j) = mean(data.spikeSnippets(:,nonArtifactIdx2),2);
                artifactAvg2(:,j)    = mean(data.spikeSnippets(:,artifactIdx2),2);
                
                clear data
                
                toc
                
            end
            
        else
            chanDataRef = load(['/media/nasko/WD_BLACK3/RatU_temp/rawDataCh' num2str(ch1) '.mat'],'data');

            spikeTimeIndxCell1 = round(SyncSp(:,1)*fs) + 1;
            [artifactAvg1,artifactWave1] = waveformAvg(double(chanDataRef.data),spikeTimeIndxCell1,36,36,300,30000,false);

            spikeTimeIndxCell1 = round(res1*fs) + 1;
            [waveAvg1,wave1] = waveformAvg(double(chanDataRef.data),spikeTimeIndxCell1,36,36,300,30000,false);

            chanDataRef = load(['/media/nasko/WD_BLACK3/RatU_temp/rawDataCh' num2str(ch2) '.mat'],'data');

            spikeTimeIndxCell2 = round(SyncSp(:,2)*fs) + 1;
            [artifactAvg2,artifactWave2] = waveformAvg(double(chanDataRef.data),spikeTimeIndxCell2,36,36,300,30000,false);

            spikeTimeIndxCell2 = round(res2*fs) + 1;
            [waveAvg2,wave2] = waveformAvg(double(chanDataRef.data),spikeTimeIndxCell2,36,36,300,30000,false);

            deltaArtifact = artifactWave1-artifactWave2;
            deltaWave     = waveAvg1 - waveAvg2;

            corrMat = corrcoef([deltaArtifact deltaWave]).^2;
        end
        
    end
    
end