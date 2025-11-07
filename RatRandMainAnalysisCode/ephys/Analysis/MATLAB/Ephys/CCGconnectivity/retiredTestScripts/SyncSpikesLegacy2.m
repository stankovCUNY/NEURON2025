function [InfoSp,InfoCCG] = SyncSpikes(n1,n2,LOI,binSize,lag)
    
    % All inputs are assumped in seconds except for lag and bin size which
    % are interger lengths of samples.
    % INPUT:
    % n1 - neuron 1 spike times
    % n2 - neuron 2 spike times
    % LOI - lags of interest, AKA significant bins 
    % binary: 1 - significant, 0 - ignore
    % lag: labels for lags, right now only using the lag size
    %
    % OUTPUT:
    % InfoSp  - informed spikes, AKA synchronous spikes
    % InfoCCG - bin counts only for significant bins. I use this for sanity
    % checks to compare with full CCG.
    
    % Atanas Stankov CUNY 2021  
    
    partition = 50; % make this adaptive in the future
    
    InfoCCG = zeros(length(lag),1);
    ISP{length(lag)} = []; % informed spike pairs per lag
    ISP = ISP';
    
    %     lagLoop = lag; % comment in to produce complete CCG
    %     lagSig = 1:length(lag);
    lagLoop = lag(LOI); % comment out to produce complete CCG
    lagSig = find(LOI); % comment out to produce complete CCG
    
    if isempty(lagLoop)
        
        InfoSp = []; % informed spikes
        return;
    else

        fs = 30000;
        n1transform  = round(n1*fs);
        n2transform  = round(n2*fs);
        
        n2LagMatrix      = repmat(n2         ,1,length(lagSig));
        n2LagMatrixTrans = repmat(n2transform,1,length(lagSig));
        
        lowerBound = lagLoop - binSize/2;
        upperBound = lagLoop + binSize/2;
        
        lowerBound = repmat(lowerBound,length(n2transform),1);
        upperBound = repmat(upperBound,length(n2transform),1);
        
%         counts = zeros(1,length(lagLoop));
        
        InfSpTemp_n1{length(lagLoop)} = [];
        InfSpTemp_n2{length(lagLoop)} = [];
        
        tic
        for i = 1:ceil(length(n1)/partition)
            
            % reset these matrices because there is recursion
            n1_3D         = [];
            n2_3D         = [];
            lowerBound_3D = [];
            upperBound_3D = []; 
            
            if i == ceil(length(n1)/partition) % last partition
                n1_partition = n1transform(1+(i-1)*partition:end);
                n1_3D = repmat(n1_partition,1,length(n2),length(lagLoop));
                n1_3D = permute(n1_3D,[2 3 1]);
                
                n2_3D      = repmat(n2LagMatrixTrans,1,1,length(n1) - (1+(i-1)*partition));
                lowerBound_3D = repmat(lowerBound,      1,1,length(n1) - (1+(i-1)*partition));
                upperBound_3D = repmat(upperBound,      1,1,length(n1) - (1+(i-1)*partition));
            else
                n1_partition = n1transform(1+(i-1)*partition:i*partition);
                n1_3D = repmat(n1_partition,1,length(n2),length(lagLoop));
                n1_3D = permute(n1_3D,[2 3 1]);

                n2_3D      = repmat(n2LagMatrixTrans,1,1,partition);
                lowerBound_3D = repmat(lowerBound,      1,1,partition);
                upperBound_3D = repmat(upperBound,      1,1,partition);
            end
            
            tempIdx    = (lowerBound_3D <= (n2_3D - n1_3D)) & ...
                         ((n2_3D - n1_3D) < upperBound_3D);
            
%             if sum(tempIdx(:)) > 0
%                 
%                 [RowIdx,ColIdx] = find(tempIdx);
%                                 
%                 for j = 1:length(ColIdx)
%                     InfSpTemp_n1{ColIdx(j)} = [InfSpTemp_n1{ColIdx(j)}; n1(i)];
%                     InfSpTemp_n2{ColIdx(j)} = [InfSpTemp_n2{ColIdx(j)}; n2LagMatrix(RowIdx(j),ColIdx(j))];
%                 end
%                 
%             end
        end
        toc
        
        for i = 1:length(lagLoop)
            ISPoi{i} = [InfSpTemp_n1{i} InfSpTemp_n2{i}];
        end
         

        InfoSp = []; % informed spikes

        for i = 1:length(ISPoi)
            InfoSp = [InfoSp; ISPoi{i}];
            ISP{lagSig(i)} = ISPoi;
            InfoCCG(lagSig(i)) = size(ISPoi{i},1);
        end

        InfoSp = min(InfoSp,[],2);
        InfoSp = unique(InfoSp);
        
    end
end