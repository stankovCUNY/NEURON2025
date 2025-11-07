function BapunCurateNoiseClusters(RatID)

    fs = 30000;
    dataDir = ['/media/nasko/3ce6d8b6-2570-4c8f-941f-9b4b0641301c/home/nasko/CUNY_Work_Han_Kamran_Utku_dataV2/data/Rat' RatID '/'];
    
    %%

    if strcmp(RatID,'N')

        neurons       = load('/media/nasko/WD_BLACK/BapunRawMultiArray/RatN/RatN_Day2_2019-10-11_03-58-54.neurons.mat');
        
        for i = 1:8
            
            spikeTimesRaw = double(readNPY(['/media/nasko/WD_BLACK/BapunRawMultiArray/RatN/RatNDay2NSD/shank' num2str(i) '/spike_times.npy']))/fs;

            spikeTimesVetted = [];
            spikeTimesMUAs   = [];

            for j = 1:length(neurons.neuron_ids) 
                if (double(neurons.shank_ids(j)) + 1) == i
                    spikeTimesVetted = [neurons.spiketrains{1,j}'; spikeTimesVetted];

                    if strcmp(neurons.neuron_type(j,:),'mua  ')
                        spikeTimesMUAs = [neurons.spiketrains{1,j}'; spikeTimesMUAs];
                    end
                end
            end

            %%

            spikeTimesNoise = setdiff(spikeTimesRaw,spikeTimesVetted);

            save([dataDir 'Bapuns_Rat' RatID '_shank_' num2str(i) '_spikeTimesNoise.mat'],'spikeTimesNoise');
            save([dataDir 'Bapuns_Rat' RatID '_shank_' num2str(i) '_spikeTimesMUAs.mat'],'spikeTimesMUAs');
        end
        
    elseif strcmp(RatID,'S')

        neurons       = load('/media/nasko/WD_BLACK/BapunRawMultiArray/RatS/RatS-Day2NSD-2020-11-27_10-22-29.neurons.mat');
        spikeTimesRaw = double(readNPY('/media/nasko/WD_BLACK/BapunRawMultiArray/RatS/RatSDay2NSD/RatS-Day2NSD-2020-11-27_10-22-29-1.GUI/spike_times.npy'))/fs;
        
        spikeTimesVetted = []; 
        spikeTimesMUAs   = [];

        for i = 1:length(neurons.neuron_ids) 
            spikeTimesVetted = [neurons.spiketrains{1,i}'; spikeTimesVetted];

            if strcmp(neurons.neuron_type(i,:),'mua  ')
                spikeTimesMUAs = [neurons.spiketrains{1,i}'; spikeTimesMUAs];
            end
        end

        %%

        spikeTimesNoise = setdiff(spikeTimesRaw,spikeTimesVetted);

        save([dataDir 'Bapuns_Rat' RatID '_spikeTimesNoise.mat'],'spikeTimesNoise');
        save([dataDir 'Bapuns_Rat' RatID '_spikeTimesMUAs.mat'],'spikeTimesMUAs');
        
    elseif strcmp(RatID,'U')

        neurons       = load('/media/nasko/WD_BLACK/BapunRawMultiArray/RatU/RatU_Day2NSD_2021-07-24_08-16-38.neurons.mat');
        spikeTimesRaw = double(readNPY('/media/nasko/WD_BLACK/BapunRawMultiArray/RatU/RatUDay2NSD/RatU_Day2NSD_2021-07-24_08-16-38-1.GUI/spike_times.npy'))/fs;
        
        spikeTimesVetted = []; 
        spikeTimesMUAs   = [];

        for i = 1:length(neurons.neuron_ids) 
            spikeTimesVetted = [neurons.spiketrains{1,i}'; spikeTimesVetted];

            if strcmp(neurons.neuron_type(i,:),'mua  ')
                spikeTimesMUAs = [neurons.spiketrains{1,i}'; spikeTimesMUAs];
            end
        end

        %%

        spikeTimesNoise = setdiff(spikeTimesRaw,spikeTimesVetted);

        save([dataDir 'Bapuns_Rat' RatID '_spikeTimesNoise.mat'],'spikeTimesNoise');
        save([dataDir 'Bapuns_Rat' RatID '_spikeTimesMUAs.mat'],'spikeTimesMUAs');
        
    end

    
    
end
