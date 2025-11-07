function  highpassFilterUtkuData(dayNo)
    
    UtkuPathRaw      = '/media/nasko/WD_BLACK/UtkuRawMultiArray/AG/';
    UtkuPathFiltered = '/media/nasko/WD_BLACK2/UtkuFilteredDataPerChannel/AG/';
    
    if str2double(dayNo) == 1
        filepath = [UtkuPathRaw 'AG_2019-12_23_NSD/'];
        savepath = [UtkuPathFiltered 'AG_2019-12_23_NSD/'];
    elseif str2double(dayNo) == 2
        filepath = [UtkuPathRaw 'AG_2019-12_27_NSD/'];
        savepath = [UtkuPathFiltered 'AG_2019-12_27_NSD/'];
    end  
    
    shankLabel    = [ones(32,1); 2*ones(32,1); 3*ones(32,1); 4*ones(32,1); 5*ones(32,1); 6*ones(32,1)];
    channelLayout = {[1:32],[33:64],[65:96],[97:128],[129:160],[161:192]}; 
    
    % generate high-pass filter
    Fstop = 250;
    Fpass = 300;
    Astop = 65;
    Apass = 0.05;
    fs    = 30000;
    d = designfilt('highpassfir','StopbandFrequency',Fstop, ...
    'PassbandFrequency',Fpass,'StopbandAttenuation',Astop, ...
    'PassbandRipple',Apass,'SampleRate',fs,'DesignMethod','equiripple');

    for channel = 1:192
        
        tic
        
        load([filepath 'shank' num2str(shankLabel(channel)) '/' 'rawDataCh' num2str(channel) '.mat'],'data'); % channel 
        data = double(data);
        
        data = filter(d,data);
        data = int16(data);  % changle from double to int16 precision
        data = data(940:end); % adjust 940 sample phase delay
        
        data = [data; zeros(940,1)];        
        
        saveName = ['ch' num2str(channel) 'highpass' num2str(Fpass) 'hz.mat'];
        save([savepath 'shank' num2str(shankLabel(channel)) '/' saveName],'data')
        
        clear data
        
        toc
        
    end

end
