fs         = 3e4;
BinSize    = 1/3e4;
Duration   = 0.002;
lagLimit   = Duration/2;

UL = 3600*fs;
LL = 0;
numElements = UL-LL+1;

r = randperm(numElements, 10);

res1 = sort(randperm(numElements, 1e6)/fs)';
res2 = sort(randperm(numElements, 1e6)/fs)';
% 
% res1 = sort(randi([0 3600*fs],1e6,1)/fs);
% res2 = sort(randi([0 3600*fs],1e5,1)/fs);

%%
tic
[ccg] = CCG([res1;res2],[ones(size(res1));2*ones(size(res2))], ...
'binSize', BinSize, 'duration', Duration, 'Fs', 1/fs, ...
'norm', 'counts');
toc

%%
tic
[CCGtracked,spikeTimesInBin] = CCGtrackedSpikeTimes(res1,res2,BinSize,lagLimit);
toc

%%
plot(ccg(:,1,2))
hold on
plot(CCGtracked)
hold off