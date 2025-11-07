function bandpassFilterHiroData(freqRange)

%     if strcmp(session_name,'RoyMaze1')
%         filepath = '/media/nasko/3ce6d8b6-2570-4c8f-941f-9b4b0641301c/home/nasko/CUNY_Work_Han_Kamran_Hiro_Data/Roy_ncs_data/maze/';
%     elseif strcmp(session_name,'KevinMaze1')
%         filepath = '/media/nasko/WD_BLACK/Kevin_ncs_data/maze/';
%     end
    
%     fpass = 300;
    fs    = 30000;

    % generate band-pass filter
    if strcmp(freqRange,'theta')
       fLow  = 5;
       fHigh = 12;
       
       savepath = '/media/nasko/WD_BLACK3/HiroRawDataPerChannelFiltered/maze/bandpassTheta/';
       
    elseif strcmp(freqRange,'gamma')
       fLow  = 25;
       fHigh = 90;
       
       savepath = '/media/nasko/WD_BLACK3/HiroRawDataPerChannelFiltered/maze/bandpassGamma/';
       
    elseif strcmp(freqRange,'ripple')
       fLow  = 130;
       fHigh = 230;
       
       savepath = '/media/nasko/WD_BLACK3/HiroRawDataPerChannelFiltered/maze/bandpassRipple/';
       
    elseif strcmp(freqRange,'agmon')
       fLow  = 350;
       fHigh = 450;
       
       savepath = '/media/nasko/WD_BLACK3/HiroRawDataPerChannelFiltered/maze/bandpassAgmon/';
       
    end
        
    fNrmLow   = fLow  / (fs/2);
    fNrmHigh  = fHigh / (fs/2);

    % determine filter coefficients:
    [b, a, k] = butter(9,[fNrmLow fNrmHigh],'bandpass');
    [sos_band,g]  = zp2sos(b, a, k);
    
    filepath = '/media/nasko/3ce6d8b6-2570-4c8f-941f-9b4b0641301c/home/nasko/CUNY_Work_Han_Kamran_Hiro_Data/Roy_ncs_data/maze/';
    
    for channel = 1:64

        [ncs] = read_neuralynx_ncs([filepath 'CSC' num2str(channel) '.ncs']); % channel 

        chanData = ncs.dat; % 512 samples per one time period which spans (1/30) * 512 ms
        onsetTime = double(ncs.TimeStamp(1));

        clear ncs

        chanData = reshape(chanData,[],1);
        chanData = filtfilt(sos_band,g,chanData);

        data.channel   = chanData;
        data.onsetTime = onsetTime;

        saveName = ['ch' num2str(channel) 'bandpass' freqRange 'range.mat'];
        save([savepath saveName],'data')

        clear chanData
        clear data

    end

end