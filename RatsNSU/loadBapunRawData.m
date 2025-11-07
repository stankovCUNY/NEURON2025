start = 0;

%% rat N

neuronsN    = load('/media/nasko/WD_BLACK/BapunRawMultiArray/RatN/RatN_Day2_2019-10-11_03-58-54.neurons.mat');
paradigmN   = load('/media/nasko/WD_BLACK/BapunRawMultiArray/RatN/RatN_Day2_2019-10-11_03-58-54.paradigm.mat');
probegroupN = load('/media/nasko/WD_BLACK/BapunRawMultiArray/RatN/RatN_Day2_2019-10-11_03-58-54.probegroup.mat');

samples = double(paradigmN.stop(3)*neuronsN.sampling_rate);
channel = double(cell2mat(probegroupN.channel_id));

rawDataPath = '/media/nasko/WD_BLACK/BapunRawMultiArray/RatN/RatN-Day2-2019-10-11_03-58-54_includes_noisy_timepoints.dat';
outDataPath = '/media/nasko/WD_BLACK2/BapunRawDataPerChannel/RatN';

% loop channels
for j = 1:length(channel)
    tic

    data = bz_LoadBinary(rawDataPath, ...
                          'samples',      samples, ...
                          'frequency',    30000, ...
                          'nChannels',    134, ...
                          'channels' ,    channel(j)); % NAT: The natural channel order should give you the channel.  E.g if spyking circus or phy files or the npy file says the max channel is 23, you would grab channel 23 from the dat file.

    saveName = ['rawDataCh' num2str(channel(j)) '.mat'];
    save([outDataPath '/' saveName],'data')

    toc
end

%% rat S

neuronsS    = load('/media/nasko/WD_BLACK/BapunRawMultiArray/RatS/RatS-Day2NSD-2020-11-27_10-22-29.neurons.mat');
paradigmS   = load('/media/nasko/WD_BLACK/BapunRawMultiArray/RatS/RatS-Day2NSD-2020-11-27_10-22-29.paradigm.mat');
probegroupS = load('/media/nasko/WD_BLACK/BapunRawMultiArray/RatS/RatS-Day2NSD-2020-11-27_10-22-29.probegroup.mat');

samples = double(paradigmS.stop(4)*neuronsS.sampling_rate);
channel = double(cell2mat(probegroupS.channel_id));

rawDataPath = '/media/nasko/WD_BLACK/BapunRawMultiArray/RatS/RatS-Day2NSD-2020-11-27_10-22-29.dat';
outDataPath = '/media/nasko/WD_BLACK2/BapunRawDataPerChannel/RatS';

% loop channels
for j = 1:length(channel)
    tic

    data = bz_LoadBinary(rawDataPath, ...
                          'samples',      samples, ...
                          'frequency',    30000, ...
                          'nChannels',    195, ...
                          'channels' ,    channel(j)); % NAT: The natural channel order should give you the channel.  E.g if spyking circus or phy files or the npy file says the max channel is 23, you would grab channel 23 from the dat file.

    saveName = ['rawDataCh' num2str(channel(j)) '.mat'];
    save([outDataPath '/' saveName],'data')

    toc
end

%% rat U

neuronsU    = load('/media/nasko/WD_BLACK/BapunRawMultiArray/RatU/RatU_Day2NSD_2021-07-24_08-16-38.neurons.mat');
paradigmU   = load('/media/nasko/WD_BLACK/BapunRawMultiArray/RatU/RatU_Day2NSD_2021-07-24_08-16-38.paradigm.mat');
probegroupU = load('/media/nasko/WD_BLACK/BapunRawMultiArray/RatU/RatU_Day2NSD_2021-07-24_08-16-38.probegroup.mat');

samples = double(paradigmU.stop(4)*neuronsU.sampling_rate);
channel = double(cell2mat(probegroupU.channel_id));

rawDataPath = '/media/nasko/WD_BLACK/BapunRawMultiArray/RatU/RatU_Day2NSD_2021-07-24_08-16-38.dat';
outDataPath = '/media/nasko/WD_BLACK2/BapunRawDataPerChannel/RatU';

% loop channels
for j = 1:length(channel)
    tic

    data = bz_LoadBinary(rawDataPath, ...
                          'samples',      samples, ...
                          'frequency',    30000, ...
                          'nChannels',    192, ...
                          'channels' ,    channel(j)); % NAT: The natural channel order should give you the channel.  E.g if spyking circus or phy files or the npy file says the max channel is 23, you would grab channel 23 from the dat file.

    saveName = ['rawDataCh' num2str(channel(j)) '.mat'];
    save([outDataPath '/' saveName],'data')

    toc
end