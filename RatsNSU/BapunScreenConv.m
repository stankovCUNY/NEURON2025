function BapunScreenConv(RatID,nullFlag)
    
    data_dir = ['/media/nasko/3ce6d8b6-2570-4c8f-941f-9b4b0641301c/home/nasko/CUNY_Work_Han_Kamran_Utku_dataV2/data/Rat' RatID '/'];

    if strcmp(RatID,'N')
        
        neurons    = load('/media/nasko/WD_BLACK/BapunRawMultiArray/RatN/RatN_Day2_2019-10-11_03-58-54.neurons.mat');
        probegroup = load('/media/nasko/WD_BLACK/BapunRawMultiArray/RatN/RatN_Day2_2019-10-11_03-58-54.probegroup.mat');

    elseif strcmp(RatID,'S')
        
        neurons    = load('/media/nasko/WD_BLACK/BapunRawMultiArray/RatS/RatS-Day2NSD-2020-11-27_10-22-29.neurons.mat');
        probegroup = load('/media/nasko/WD_BLACK/BapunRawMultiArray/RatS/RatS-Day2NSD-2020-11-27_10-22-29.probegroup.mat');

    elseif strcmp(RatID,'U')
        
        neurons    = load('/media/nasko/WD_BLACK/BapunRawMultiArray/RatU/RatU_Day2NSD_2021-07-24_08-16-38.neurons.mat');
        probegroup = load('/media/nasko/WD_BLACK/BapunRawMultiArray/RatU/RatU_Day2NSD_2021-07-24_08-16-38.probegroup.mat');
        
    end
    
    spikeTimes = []; % spiket for Eran
    spikeInd   = []; % spikeind for Eran
    
    for i = 1:length(neurons.neuron_ids) 
        spikeTimes = [neurons.spiketrains{1,i}'; spikeTimes];
        spikeInd   = [ones(length(neurons.spiketrains{1,i}),1)*(double(neurons.neuron_ids(i)) + 1); spikeInd];
    end
    
    if nullFlag
        one_ms = 0.001;
        delta = 5*one_ms;
        spikeTimes = floor(spikeTimes/delta)*delta + rand(size(spikeTimes))*delta;
    end
    
    cellID     = double(neurons.neuron_ids) + 1;    % Cells for Eran
    shankID    = double(neurons.shank_ids) + 1;     % shank for Eran
    SampleRate = 30000;                             % SampleRate for Eran
    alpha      = 0.05;                              % alpah for Eran
    wintype    = 'gauss';                           % wintype for Eran
    
    %% jscale 1
    jscale = 1; % jscale for Eran

    [ExcPairs, ...
     InhPairs, ...
     GapPairs, ...
     Rzero] = EranConv_group_forUtku(spikeTimes, ...
                                     spikeInd, ...
                                     cellID, ...
                                     SampleRate, ...
                                     jscale, ...
                                     alpha,...
                                     shankID, ...
                                     wintype);

    pairs.ExcPairs   = ExcPairs;
    pairs.InhPairs   = InhPairs;
    pairs.GapPairs   = GapPairs;
    pairs.RZero      = Rzero;
    pairs.jscale     = jscale;
    pairs.nullFlag   = nullFlag;
    pairs.spikeTimes = spikeTimes;
    pairs.spikeInd   = spikeInd;
    
    if nullFlag
        save(fullfile(data_dir,['Bapuns_Rat' RatID 'null_jscale' num2str(jscale) '_alpha' ...
                    num2str(round(alpha*100)) '_pairs']), ...
                    'pairs', 'jscale', 'alpha')
    else 
        save(fullfile(data_dir,['Bapuns_Rat' RatID '_jscale' num2str(jscale) '_alpha' ...
                    num2str(round(alpha*100)) '_pairs']), ...
                    'pairs', 'jscale', 'alpha')
    end
    
    %% jscale 5
    jscale = 5;  % jscale for Eran

    [ExcPairs, ...
     InhPairs, ...
     GapPairs, ...
     Rzero] = EranConv_group_forUtku(spikeTimes, ...
                                     spikeInd, ...
                                     cellID, ...
                                     SampleRate, ...
                                     jscale, ...
                                     alpha,...
                                     shankID, ...
                                     wintype);

    pairs.ExcPairs   = ExcPairs;
    pairs.InhPairs   = InhPairs;
    pairs.GapPairs   = GapPairs;
    pairs.RZero      = Rzero;
    pairs.jscale     = jscale;
    pairs.nullFlag   = nullFlag;
    pairs.spikeTimes = spikeTimes;
    pairs.spikeInd   = spikeInd;
    
    if nullFlag
        save(fullfile(data_dir,['Bapuns_Rat' RatID 'null_jscale' num2str(jscale) '_alpha' ...
                    num2str(round(alpha*100)) '_pairs']), ...
                    'pairs', 'jscale', 'alpha')
    else
        save(fullfile(data_dir,['Bapuns_Rat' RatID '_jscale' num2str(jscale) '_alpha' ...
                    num2str(round(alpha*100)) '_pairs']), ...
                    'pairs', 'jscale', 'alpha')
    end
 end