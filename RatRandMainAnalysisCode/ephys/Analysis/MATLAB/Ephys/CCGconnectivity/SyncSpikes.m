function [SyncSp,SyncCCG,SyncSpBin,SyncSpIdxBin] = SyncSpikes(n1,n2,LOI,duration)

    % Run the CCG function:
    [tsOffsets, n1idx, n2idx] = crosscorrelogram(n1, n2, [-(duration/2),duration/2]);

    % To make a CCG:
%     lgs = -105*(1/30000):1/30000:105*(1/30000);
%     CCG = histcounts(tsOffsets*30000,'NumBins',211,'BinMethod','integers');

    n1entry = n1(n1idx);
    n2entry = n2(n2idx);
    
    n1SpkIdx = 1:length(n1);
    n2SpkIdx = 1:length(n2);
    
    n1SpkIdx = n1SpkIdx(n1idx)';
    n2SpkIdx = n2SpkIdx(n2idx)';
    
    SyncSpBin = {};
    
    tsOffsets = round(tsOffsets*30000);
    lagInt = -((length(LOI)-1)/2):((length(LOI)-1)/2);
    
    for i = 1:length(LOI)
        idx = find(tsOffsets == lagInt(i));
        SyncSpBin{i}    = [n1entry(idx)  n2entry(idx)];
        SyncSpIdxBin{i} = [n1SpkIdx(idx) n2SpkIdx(idx)];
    end
    
    SyncCCG = zeros(length(LOI),1);
    LOIidx = find(LOI);
    SyncSp = [];
    
    for i = 1:length(LOIidx)
        SyncSp = [SyncSp; SyncSpBin{LOIidx(i)}];
        SyncCCG(LOIidx(i)) = size(SyncSpBin{LOIidx(i)},1);
    end
    
%     SyncSp = min(SyncSp,[],2);
end