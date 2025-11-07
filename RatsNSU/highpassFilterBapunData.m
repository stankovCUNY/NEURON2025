function  highpassFilterBapunData
    
    animalIDs  = {'N','S','U'};
    numOfChans = [128,195,192];
    
    for loopAnimal = 1:length(animalIDs)
        
        filepath = ['/media/nasko/WD_BLACK21/BapunRawDataPerChannel/Rat'     animalIDs{loopAnimal} '/'];
        savepath = ['/media/nasko/WD_BLACK3/BapunFilteredDataPerChannel/Rat' animalIDs{loopAnimal} '/'];

        % generate high-pass filter
        Fstop = 250;
        Fpass = 300;
        Astop = 65;
        Apass = 0.05;
        fs    = 30000;
        d = designfilt('highpassfir','StopbandFrequency',Fstop, ...
        'PassbandFrequency',Fpass,'StopbandAttenuation',Astop, ...
        'PassbandRipple',Apass,'SampleRate',fs,'DesignMethod','equiripple');
        
        loopLenght = numOfChans(loopAnimal);
    
        % filter channels one by one
        for loopChannel = 1:loopLenght
            
            if strcmp(animalIDs{loopAnimal},'S')
                if (loopChannel == 52) || (loopChannel >= 193)
                    continue
                end
            end
            
            tic

            load([filepath '/' 'rawDataCh' num2str(loopChannel) '.mat'],'data'); % channel 
            data = double(data);

            data = filter(d,data);
            data = int16(data);  % changle from double to int16 precision
            data = data(940:end); % adjust 940 sample phase delay

            data = [data; zeros(940,1)];        

            saveName = ['ch' num2str(loopChannel) 'highpass' num2str(Fpass) 'hz.mat'];
            save([savepath saveName],'data')

            clear data

            toc

        end
    end
end
