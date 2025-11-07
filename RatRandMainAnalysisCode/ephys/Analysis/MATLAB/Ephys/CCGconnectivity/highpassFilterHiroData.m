function  highpassFilterHiroData(fpass,session_name)
    
    if strcmp(session_name,'RoyMaze1')
        filepath = '/media/nasko/3ce6d8b6-2570-4c8f-941f-9b4b0641301c/home/nasko/CUNY_Work_Han_Kamran_Hiro_Data/Roy_ncs_data/maze/';
        fs       = 30000;
    elseif strcmp(session_name,'KevinMaze1')
        filepath = '/media/nasko/WD_BLACK3/Kevin_ncs_data/maze/';
        fs       = 32000;
    end
    
%     fpass = 300;

     % generate high-pass filter
    [b, a, k] = butter(9, fpass / fs * 2, 'high');
    [sos_high,g] = zp2sos(b, a, k);

    for channel = 1:64
        
        [ncs] = read_neuralynx_ncs([filepath 'CSC' num2str(channel) '.ncs']); % channel 

        chanData = ncs.dat; % 512 samples per one time period which spans (1/30) * 512 ms
        onsetTime = double(ncs.TimeStamp(1));
        
        clear ncs
        
        chanData = reshape(chanData,[],1);
        chanData = filtfilt(sos_high,g,chanData);
  
        data.channel   = chanData;
        data.onsetTime = onsetTime;
        
        saveName = ['ch' num2str(channel) 'highpass' num2str(fpass) 'hz.mat'];
        save([filepath saveName],'data')
        
        clear chanData
        clear data
        
    end

end