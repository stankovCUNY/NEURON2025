function BapunWaveformSnippets(RatID)   

%% spike times

    randomFlag = true; 
    
    
    if strcmp(RatID,'N')
        
        neurons      = load('/media/nasko/WD_BLACK/BapunRawMultiArray/RatN/RatN_Day2_2019-10-11_03-58-54.neurons.mat');
        rawDataPath  = '/media/nasko/WD_BLACK/BapunRawMultiArray/RatN/RatN-Day2-2019-10-11_03-58-54_includes_noisy_timepoints.dat';
        outDatapath  = '/media/nasko/WD_BLACK2/BapunRawDataPerNeuron/RatN/';
        
        % hard coding from probe file
        channelSet = {1:16, ...
                      17:32, ... 
                      33:48, ...
                      49:64, ...
                      65:80, ...
                      81:96, ...
                      97:112, ...
                      113:128};
        
        % hard coding from .xml
        NoOfchannels = 134; 
        
    elseif strcmp(RatID,'S')
        
        neurons      = load('/media/nasko/WD_BLACK/BapunRawMultiArray/RatS/RatS-Day2NSD-2020-11-27_10-22-29.neurons.mat');
        rawDataPath  = '/media/nasko/WD_BLACK/BapunRawMultiArray/RatS/RatS-Day2NSD-2020-11-27_10-22-29.dat';
        outDatapath  = '/media/nasko/WD_BLACK2/BapunRawDataPerNeuron/RatS/';
        
        % hard coding from probe file
        channelSet = {1:10, ...
                      11:21, ...
                      22:32, ... 
                      33:43, ...
                      44:54, ...
                      55:63, ...
                      64:79, ...
                      80:95, ...
                      96:111, ...
                      112:127, ...
                      128:143, ...
                      144:159, ...
                      160:175, ...
                      176:191};
        
        % hard coding from .xml
        NoOfchannels = 195; 
        
    elseif strcmp(RatID,'U')
        
        neurons      = load('/media/nasko/WD_BLACK/BapunRawMultiArray/RatU/RatU_Day2NSD_2021-07-24_08-16-38.neurons.mat');
        rawDataPath  = '/media/nasko/WD_BLACK/BapunRawMultiArray/RatU/RatU_Day2NSD_2021-07-24_08-16-38.dat';
        outDatapath  = '/media/nasko/WD_BLACK2/BapunRawDataPerNeuron/RatU/';
        
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
        
        % hard coding from .xml
        NoOfchannels = 192; 
        
    end
    
    halfWidth         = 27;
    numberOfWaveforms = 1000;
    
    % loop neurons
    for unitIdx = 1:length(neurons.neuron_ids)

        fs = double(neurons.sampling_rate);
          
        cellType   = neurons.neuron_type(unitIdx,:);
        spikeTimes = neurons.spiketrains{1,unitIdx}; % in seconds
        spikeTimes = spikeTimes*fs;                  % convert to samples
        
        % 1000 random spike times
        if randomFlag
           if numberOfWaveforms < length(spikeTimes)
                spikeTimes = spikeTimes(randperm(length(spikeTimes),numberOfWaveforms));
           end 
        end
        spikeTimes = sort(spikeTimes);
        
%         if strcmp(RatID,'N') % noisy time adjustment for rat N
%             noisyIdx = find(spikeTimes > (168*60*fs));               % units: [min * sec * smaples]
%             spikeTimes(noisyIdx) = spikeTimes(noisyIdx) + (4*60*fs); % units: [min * sec * smaples]
%         end
        
        % channel    = 1:128;
        
        % channel    = double(neurons.peak_channels(1)) + 1;

        if strcmp(cellType,'pyr  ')
            cellType = 'p';
        elseif strcmp(cellType,'inter')
            cellType = 'i';
        elseif strcmp(cellType,'mua  ')
            disp('mua unit skipped!')
            continue
        end

        %% raw data

        for channelSetIdx = 1:size(channelSet,2)

            waveformSnippetMat = zeros(halfWidth*2,length(channelSet{channelSetIdx}),length(spikeTimes));
            
            snippet = struct();
            
            tic
            parfor i = 1:length(spikeTimes)

                offset = round(spikeTimes(i) - halfWidth);
                waveformSnippet = bz_LoadBinary(rawDataPath, ...
                                      'offset',       offset, ... 
                                      'samples',      halfWidth*2, ...
                                      'frequency',    fs, ...
                                      'nChannels',    NoOfchannels, ...
                                      'channels' ,    channelSet{channelSetIdx});

               waveformSnippetMat(:,:,i) = waveformSnippet;
            end
            toc

            if randomFlag
                fileName = ['neuron' num2str(neurons.neuron_ids(unitIdx) + 1) cellType 'Shank' num2str(channelSetIdx) 'random' num2str(numberOfWaveforms) 'spikes'];
            else
                fileName = ['neuron' num2str(neurons.neuron_ids(unitIdx) + 1) cellType 'Shank' num2str(channelSetIdx)];
            end

            snippet.waveformSnippetMat = waveformSnippetMat;
            snippet.spikeTimes         = spikeTimes;

            save([outDatapath fileName],'snippet')
        end

        disp([num2str(100*round(unitIdx/length(neurons.neuron_ids),4)) '% done'])

    end 
end