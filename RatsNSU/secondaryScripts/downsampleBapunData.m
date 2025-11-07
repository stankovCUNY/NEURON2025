function downsampleBapunData
    
    filepath = '/media/nasko/WD_BLACK3/RatU_temp/highpassFiltered/';
    savepath = '/media/nasko/WD_BLACK3/RatU_temp/highpassFiltered/downsampled1250hz/';

    for channel = 1:192
        
        tic
        
        load([filepath 'ch' num2str(channel) 'highpass300hz.mat'],'data'); % channel 
        data = double(data);
        
        ratio = 30000/1250; % down sample from 30khz to 1250hz
        
        data = downsample(data,ratio);

        saveName = ['ch' num2str(channel) 'downsampled1250hz.mat'];
        save([savepath saveName],'data')
        
        clear data
        
        toc
        
    end

end
