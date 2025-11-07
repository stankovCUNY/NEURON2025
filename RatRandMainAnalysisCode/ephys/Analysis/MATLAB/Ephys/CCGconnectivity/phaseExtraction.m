function phase = phaseExtraction(ch,freqRange,session_name,fs)
    
%     session_name   = 'RoyMaze1';
    
    if strcmp(freqRange,'theta')
       fLow  = 5;
       fHigh = 12;
       
       loadpath = '/media/nasko/WD_BLACK3/HiroRawDataPerChannelFiltered/maze/bandpassTheta/';
       
    elseif strcmp(freqRange,'gamma')
       fLow  = 25;
       fHigh = 90;
       
       loadpath = '/media/nasko/WD_BLACK3/HiroRawDataPerChannelFiltered/maze/bandpassGamma/';
       
    elseif strcmp(freqRange,'ripple')
       fLow  = 130;
       fHigh = 230;
       
       loadpath = '/media/nasko/WD_BLACK3/HiroRawDataPerChannelFiltered/maze/bandpassRipple/';
       
    elseif strcmp(freqRange,'agmon')
       fLow  = 350;
       fHigh = 450;
       
       loadpath = '/media/nasko/WD_BLACK3/HiroRawDataPerChannelFiltered/maze/bandpassAgmon/';
       
    end
%         
%     fNrmLow   = fLow  / (fs/2);
%     fNrmHigh  = fHigh / (fs/2);
% 
%     x = HiroLoadRawNSC(session_name,ch);
% 
%     % determine filter coefficients:
%     [b, a, k] = butter(9,[fNrmLow fNrmHigh],'bandpass');
%     [sos_band,g]  = zp2sos(b, a, k);

%     y = filtfilt(sos_band,g,x);
    
    load([loadpath 'ch' num2str(ch) 'bandpass' freqRange 'range.mat'])
    y = data.channel;
    z = hilbert(y);

%     phase = atan(imag(z)./real(z)); % What is the bug here!?
    
    phase = angle(z);
    
%     plot(x(1:3e4))
%     hold on
%     plot(y(1:3e4))
%     hold off
end