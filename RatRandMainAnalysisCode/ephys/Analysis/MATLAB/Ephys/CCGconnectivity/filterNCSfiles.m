session_name = 'RoyMaze1';

fs    = 30000; 
fpass = 300;

% generate high-pass filter
[b, a, k] = butter(9, fpass / fs * 2, 'high');
sos_high = zp2sos(b, a, k);

if strcmp(session_name,'RoyMaze1')
    filepath = '/media/nasko/3ce6d8b6-2570-4c8f-941f-9b4b0641301c/home/nasko/CUNY_Work_Han_Kamran_Hiro_Data/Roy_ncs_data/maze/';
elseif strcmp(session_name,'KevinMaze1')
    filepath = '/media/nasko/WD_BLACK/Kevin_ncs_data/maze/';
end

for channel = 1:65
    
    [ncs] = read_neuralynx_ncs([filepath 'CSC' num2str(channel) '.ncs']); % channel 

    chanData = ncs.dat; % 512 samples per one time period which spans (1/30) * 512 ms
%     chanData = reshape(chanData,[],1);
    chanData = chanData(:); 
    
    tic 
    % high-pass filter the data
%     chanDataFilt = highpass(chanData,fpass,fs);
    chanDataFilt = sosfilt(sos_high, chanData);
    
    toc 
    
%     onsetTime = double(ncs.TimeStamp(1));
    
end