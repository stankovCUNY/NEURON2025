function [InfoSp,InfoCCG] = SyncSpikesLegacy(n1,n2,LOI,binSize,lag)
    
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
    
    % flip neurons if n1 has more spikes than n2
%     if length(n1) > length(n2)
%         
%         n1temp = n1; 
%         n2temp = n2;
%         
%         n2 = n1temp; 
%         n1 = n2temp;
%         
%         clear n1temp n2temp
%         
%         LOI = flipud(LOI);
%         
%         flipped = 1;
%     else 
%         flipped = 0; 
%     end
    
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
        
        for i = 1:length(n1)
                        
            tempIdx    = (lowerBound <= (n2LagMatrixTrans - n1transform(i))) & ...
                         ((n2LagMatrixTrans - n1transform(i)) < upperBound);
            
            if sum(tempIdx(:)) > 0
                
                [RowIdx,ColIdx] = find(tempIdx);
                                
                for j = 1:length(ColIdx)
                    InfSpTemp_n1{ColIdx(j)} = [InfSpTemp_n1{ColIdx(j)}; n1(i)];
                    InfSpTemp_n2{ColIdx(j)} = [InfSpTemp_n2{ColIdx(j)}; n2LagMatrix(RowIdx(j),ColIdx(j))];
                end
                
%                 break;
            end
            
        end
        
        for i = 1:length(lagLoop)
            ISPoi{i} = [InfSpTemp_n1{i} InfSpTemp_n2{i}];
        end
         
        %% parfor method
        
%         % make InfoCCG and ISP compatible with parfor
%         InfoCCGparfor = zeros(length(lagLoop),1);
%         ISPparfor{length(lagLoop)} = []; % informed spike pairs per lag
%         ISPparfor = ISPparfor';
%         
%         parfor k = 1:length(lagLoop) % loop significant bins
% 
%             count = 0;
%             InfSpTemp = [];
% 
%             lowerBound = lagLoop(k) - binSize/2;
%             upperBound = lagLoop(k) + binSize/2;
% 
%             %% legacy method - super slow nested loop for n by m pairwise substracions
% 
%     %         for i = 1:length(n1)
%     %             for j = 1:length(n2)
%     % 
%     %                 deltaTimes = n2transform(j) - n1transform(i);
%     %                 
%     %                 
%     %                 if (lowerBound <= deltaTimes) && (deltaTimes < upperBound)
%     % 
%     %                     count = count + 1;
%     %                     InfSpTemp = [[n1(i) n2(j)]; InfSpTemp];
%     %                     
%     %                 end
%     %             end
%     %         end
% 
%             %% new method
% 
%             InfSpTemp_n1 = [];
%             InfSpTemp_n2 = [];
% 
%             for i = 1:length(n1)
% 
%                 deltaTimes = n2transform - n1transform(i);
% 
%                 tempIdx = (lowerBound <= deltaTimes) & (deltaTimes < upperBound);
% 
%                 count = count + sum(tempIdx);
% 
%                 if sum(tempIdx) > 0
%                     InfSpTemp_n1 = [InfSpTemp_n1 n1(i)];
%                     InfSpTemp_n2 = [InfSpTemp_n2 n2(tempIdx)];
%                 end
% 
%                 InfSpTemp = [InfSpTemp_n1; InfSpTemp_n2]';
%             end
% 
%             ISPparfor{k} = InfSpTemp;
%             InfoCCGparfor(k) = count;
%         end

%         for k = 1:length(lagLoop) % loop significant bins
% 
%             ISP{lagSig(k)} = ISPparfor{k};
%             InfoCCG(lagSig(k)) = InfoCCGparfor(k);
% 
%         end 

%%

%         ISPoi = ISP(LOI); % informed spike pairs per lag of interest
        InfoSp = []; % informed spikes

        for i = 1:length(ISPoi)
            InfoSp = [InfoSp; ISPoi{i}];
            ISP{lagSig(i)} = ISPoi;
            InfoCCG(lagSig(i)) = size(ISPoi{i},1);
        end

        InfoSp = min(InfoSp,[],2);
        InfoSp = unique(InfoSp);
        
%         InfoCCG(lagSig) = counts;
        
%         if flipped
%             InfoCCG = flipud(InfoCCG);
%         end
    end
end