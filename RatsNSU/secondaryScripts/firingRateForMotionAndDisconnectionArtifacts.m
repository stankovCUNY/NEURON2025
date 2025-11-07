clear all

%% rat U

neuronsU    = load('/media/nasko/WD_BLACK3/RatU_temp2/RatU_Day2NSD_2021-07-24_08-16-38.neurons.mat');
paradigmU   = load('/media/nasko/WD_BLACK3/RatU_temp2/RatU_Day2NSD_2021-07-24_08-16-38.paradigm.mat');
probegroupU = load('/media/nasko/WD_BLACK3/RatU_temp2/RatU_Day2NSD_2021-07-24_08-16-38.probegroup.mat');

samples = double(paradigmU.stop(4)*neuronsU.sampling_rate);
channel = double(cell2mat(probegroupU.channel_id));

spikeBin = zeros(samples,1);
neuronCount = 0;

fs = 30000;

winSize = fs;
% winSize = 0;

for i = 1:length(neuronsU.spiketrains)
    
    i
    
    cellType = neuronsU.neuron_type(i,:);
    
    % reject MUAs at this stage
    if strcmp(cellType,'mua  ')
        continue
    else
        neuronCount = neuronCount + 1;
    end
    
    %% curate spike times

    fs   = 30000;

    res  = neuronsU.spiketrains{1,i}*fs; % change units from seconds to samples
    res  = round(res) + 1;               % correct for miniscule decimals bc of change of units
    res(res > samples) = [];
    
    %% estimating firing rates

    spikeBin(res) = spikeBin(res) + 1;
    spikeBin = single(spikeBin);
end

if winSize == 0
    spikeCount = spikeBin/neuronCount;
else
    spikeCount = movsum(spikeBin/neuronCount,winSize); % how to normalize here?
end

clear spikeBin

cutoff = mean(spikeCount) + std(spikeCount);
% cutoff = median(spikeCount);

spikeCount(spikeCount<cutoff) =  0;

filterRes = find(spikeCount == 0);