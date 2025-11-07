function [spindleMaze,rippleMaze,mobileMaze,thetaMaze,gammaMaze] = makeBXandFreqFlags
    
    dataPath = '/home/nasko/CUNY_Work_Han_Kamran_Nat/data/wake_new/';

    %% session markers 
    load([dataPath 'wake-behavior.mat'])
    mazeOn  = behavior.RoyMaze1.time(2,1); 
    mazeOff = behavior.RoyMaze1.time(2,2);

    %% ripple
    load([dataPath 'wake-ripple.mat'])
    rippleMaze    = ripple.RoyMaze1.time;
    rippleMazeBin = (rippleMaze < mazeOff) & (rippleMaze > mazeOn);
    rippleMaze    = rippleMaze(rippleMazeBin(:,1) & rippleMazeBin(:,2),:);

    %% spindle
    load([dataPath 'wake-spindle.mat'])
    spindleMaze    = spindle.RoyMaze1.HPC.time;
    spindleMazeBin = (spindleMaze < mazeOff) & (spindleMaze > mazeOn);
    spindleMaze    = spindleMaze(spindleMazeBin(:,1) & spindleMazeBin(:,2),:);

    %% speed

    load([dataPath 'wake-speed.mat'])

    v = speed.RoyMaze1.v;
    t = speed.RoyMaze1.t; % fs = 100hz, raw time data is in microsecond units
    win_size = 1000; % 10 sec guass window length, how much should I smooth??
    mobileThr = 2;

    win = gausswin(win_size); % create Gaussian kernel
    win = win/sum(win); % normalize the Gaussian kernel

    v_filt = filtfilt(win,1,v); % this function will do the filtering with no introduced lag
    v_filt_bin = zeros(length(v),1); % create vector of zeros with length of total recording session 
    v_filt_bin(find(v_filt >= mobileThr)) = 1; % threshold the filtered times series to create a binary event matrix

    onset_mobile          = t(find(diff(v_filt_bin) == 1)  + 1);  % triggers for mobility onset
    offset_mobile         = t(find(diff(v_filt_bin) == -1) + 1); % triggers for mobility offset

    mobileMaze    = [onset_mobile offset_mobile];
    mobileMazeBin = (mobileMaze < mazeOff) & (mobileMaze > mazeOn);
    mobileMaze    = mobileMaze(mobileMazeBin(:,1) & mobileMazeBin(:,2),:);

    %% spectrum

    load('wake-spectrum.mat')

    Pxx = spectrum.RoyMaze1.Pxx;
    t   = spectrum.RoyMaze1.time; % fs = 2.4414hz

    thetaPxx = mean(Pxx(:,10:21),2);
    gammaPxx = mean(Pxx(:,42:end),2);
    notPxx   = mean(Pxx(:));

    thetaLevel = 0.5*log(thetaPxx/notPxx);
    gammaLevel = 0.5*log(gammaPxx/notPxx);

    threshDB = 0.5;

    thetaLevelBin = thetaLevel > threshDB;
    gammaLevelBin = gammaLevel > threshDB;

    onset_theta  = t(find(diff(thetaLevelBin) == 1)  + 1);
    offset_theta = t(find(diff(thetaLevelBin) == -1) + 1);

    thetaMaze    = [onset_theta offset_theta];
    thetaMazeBin = (thetaMaze < mazeOff) & (thetaMaze > mazeOn);
    thetaMaze    = thetaMaze(thetaMazeBin(:,1) & thetaMazeBin(:,2),:);

    onset_gamma  = t(find(diff(gammaLevelBin) == 1)  + 1);
    offset_gamma = t(find(diff(gammaLevelBin) == -1) + 1);

    gammaMaze    = [onset_gamma offset_gamma];
    gammaMazeBin = (gammaMaze < mazeOff) & (gammaMaze > mazeOn);
    gammaMaze    = gammaMaze(gammaMazeBin(:,1) & gammaMazeBin(:,2),:);
    
    %% convert to seconds
    spindleMaze = spindleMaze/1e6;
    rippleMaze  = rippleMaze/1e6;
    mobileMaze  = mobileMaze/1e6;
    thetaMaze   = thetaMaze/1e6;
    gammaMaze   = gammaMaze/1e6;
    
end




