 function UtkuCurateNoiseClusters(dayNo,shankNo)

    fs = 30000;
    UtkuPath    = '/media/nasko/3ce6d8b6-2570-4c8f-941f-9b4b0641301c/home/nasko/CUNY_Work_Han_Kamran_Utku_dataV2/data/';
    
    day   = str2double(dayNo);
    shank = str2double(shankNo);
    
    if day == 1
        datapath = [UtkuPath 'AG_2019-12-23_NSD' '/'];
    elseif day == 2
        datapath = [UtkuPath 'AG_2019-12-27_NSD' '/'];
    end 
    
    %%

    neurons       = loadUtkuV2(dayNo);
    spikeTimesRaw = double(readNPY([datapath 'Units/s' shankNo '/spike_times.npy']))/fs;

    spikeTimesVetted = []; 

    for i = 1:length(neurons.cell_metrics.UID) 
        if neurons.cell_metrics.shankID(i) == shank
            spikeTimesVetted = [neurons.cell_metrics.spikes.times{1,i}; spikeTimesVetted];
        end
    end

    %%

    spikeTimesNoise = setdiff(spikeTimesRaw,spikeTimesVetted);
    
    save([datapath 'Utkus_RatAG_day' dayNo '_shank' shankNo '_spikeTimesNoise.mat'],'spikeTimesNoise');
    
end
