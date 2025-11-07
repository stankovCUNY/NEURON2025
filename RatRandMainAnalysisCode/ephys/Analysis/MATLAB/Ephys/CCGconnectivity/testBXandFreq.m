%% position

load('wake-position.mat')

X = position.RoyMaze1.x;
Y = position.RoyMaze1.y;

XBinEdges = [0 70:5:650 650];
YBinEdges = [160 160:5:400 400];

histogram2(X,Y,XBinEdges,YBinEdges,'Normalization','probability','FaceColor', ...
'flat','DisplayStyle','tile','ShowEmptyBins','on','LineStyle','none');

colorbar
caxis([0 0.001])

%% speed

load('wake-speed.mat')

v = speed.RoyMaze1.v;
t = speed.RoyMaze1.t; % fs = 100hz, raw time data is in microsecond units
win_size = 1000; % 10 sec guass window length, how much should I smooth??
mobileThr = 2;

win = gausswin(win_size); % create Gaussian kernel
win = win/sum(win); % normalize the Gaussian kernel

v_filt = filtfilt(win,1,v); % this function will do the filtering with no introduced lag
v_filt_bin = zeros(length(v),1); % create vector of zeros with length of total recording session 
v_filt_bin(find(v_filt >= mobileThr)) = 1; % threshold the filtered times series to create a binary event matrix

onset_mobile          = find(diff(v_filt_bin) == 1) + 1;  % triggers for mobility onset
offset_mobile         = find(diff(v_filt_bin) == -1) + 1; % triggers for mobility offset

%% spectrum

load('wake-spectrum.mat')

Pxx = spectrum.RoyMaze1.Pxx;
t   = spectrum.RoyMaze1.time; % fs = 2.4414hz

imagesc(Pxx')
set(gca, 'YTick', 1:83, 'YTickLabel', spectrum.RoyMaze1.freq)
colorbar
caxis([1e5 3e5])

thetaPxx = mean(Pxx(:,10:21),2);
gammaPxx = mean(Pxx(:,42:end),2);
notPxx   = mean(Pxx(:));

thetaLevel = 0.5*log(thetaPxx/notPxx);
gammaLevel = 0.5*log(gammaPxx/notPxx);

threshDB = 0.5;

thetaLevelBin = thetaLevel > threshDB;
gammaLevelBin = gammaLevel > threshDB;


